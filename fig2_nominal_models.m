%FIG2_NOMINAL_MODELS  Fig. 2 of Sun et al.: open-loop approximation error of the DCM-based and
%   quaternion-based nominal Koopman models, Nphi = 3, 5, 7, 9, under the Sec. V-B test torque.
%   Both models are lifted once at t = 0 and then run open loop for 300 s; the error is
%   ||C z_k - x_k||_2 against the true nonlinear dynamics (ode45).
%   Writes results/fig2_nominal_models.png.
%
%   [paper] = quoted from the paper
%   [guess] = not given in the paper, chosen here

here = fileparts(mfilename('fullpath'));
if isempty(here) || contains(here, 'Editor_'), here = pwd; end
addpath(fullfile(here, 'lib'));

% ---------------- settings ----------------
p = pakmpc_params();
nphis = [3 5 7 9];          % [paper] Fig. 2 legend
tend = 300;                 % [paper] Fig. 2 time axis
jc0 = p.jc;                 % [guess] no inertia uncertainty, as in Sec. V-B for Fig. 3
x0 = [0; 0; 0; p.q0];       % [guess] w0 = 0 matches Fig. 2's flat start; the Sec. V-A w0 gives
                            %   large errors from t = 0 (set x0 = [p.w0; p.q0] to see)
dist_on = 0;                % [guess] add tau_d of Sec. V-C to the true plant; = 1 matches the
                            %   300 s error levels but not the flat start
tau = @(t) 0.05*[7.5*sin(2*t) + 7.5*sin(0.66*t);
                 2.5*sin(2.2*t + pi/2) + 2.5*sin(0.22*t) + 1;
                 5*sin(0.66*t) + 5*sin(1.52*t) - 1];        % [paper] Sec. V-B, N m

% ---------------- truth ----------------
t = 0:p.ts:tend; nt = numel(t);
opt = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
[~, x] = ode45(@(tt, xx) attitude_ode(xx, tau(tt) + dist_on*p.dist(tt), p.jc), t, x0, opt);
x = x';

% ---------------- open-loop Koopman predictions ----------------
builders = {@dcm_koopman, @quat_koopman};
err = zeros(numel(builders), numel(nphis), nt);
for i = 1:numel(builders)
    for j = 1:numel(nphis)
        m = builders{i}(nphis(j), jc0, p.ts);
        err(i, j, :) = predict_error(m, x, t, tau);
    end
end

% ---------------- report ----------------
fprintf('%-18s %5s %9s %9s %9s\n', 'model', 'Nphi', 'e(100 s)', 'e(200 s)', 'e(300 s)');
ik = round([100 200 300]/p.ts) + 1;
for i = 1:numel(builders)
    m = builders{i}(nphis(1), jc0, p.ts);
    for j = 1:numel(nphis)
        fprintf('%-18s %5d %9.4f %9.4f %9.4f\n', m.name, nphis(j), squeeze(err(i, j, ik)));
    end
end

% ---------------- plot, styled after Fig. 2 ----------------
styles = {'-', '--', ':', '-.'};
cols = [0.85 0.2 0.1; 0.1 0.2 0.85; 0.5 0 0.5; 0.3 0.7 0.3];
titles = {sprintf('DCM-based\nKoopman model'), sprintf('Quaternion-based\nKoopman model')};
fig = figure('Position', [100 100 700 320], 'Visible', 'off');
theme(fig, 'light');
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact');
for i = 1:numel(builders)
    ax = nexttile(tl); hold(ax, 'on'); box(ax, 'on');
    for j = 1:numel(nphis)
        plot(ax, t, squeeze(err(i, j, :)), styles{j}, 'Color', cols(j, :), 'LineWidth', 1.2, ...
            'DisplayName', sprintf('N_\\phi = %d', nphis(j)));
    end
    title(ax, titles{i}, 'FontWeight', 'normal');
    xlabel(ax, 'Time (s)'); ylabel(ax, 'Approximation errors');
    xlim(ax, [0 tend]);
    legend(ax, 'Location', 'northwest');
end
out = fullfile(here, 'results');
if ~exist(out, 'dir'), mkdir(out); end
exportgraphics(fig, fullfile(out, 'fig2_nominal_models.png'), 'Resolution', 200);
fprintf('saved %s\n', fullfile(out, 'fig2_nominal_models.png'));

% ---------------- helpers ----------------
function e = predict_error(m, x, t, tau)
% Lift once, then z+ = A z + B(w_hat, q_hat) u with u = -w_hat^x Jc0 w_hat + tau (Eq. 8).
% [guess] B is re-evaluated at the model's own predicted state, so the run is a pure prediction.
dcm = strcmp(m.name, 'DCM-based');
nt = numel(t);
e = zeros(1, nt);
z = m.lift(x(:, 1));
for k = 1:nt
    if dcm
        r = q2r(x(4:7, k));
        e(k) = norm(m.c*z - [x(1:3, k); r(:)]);
        [u, ~, v] = svd(reshape(z(4:12), 3, 3));
        qh = rotm2quat(u*v')';          % nearest rotation of the predicted vec(R)
    else
        e(k) = norm(m.c*z - x(:, k));
        qh = z(4:7);
    end
    wh = z(1:3);
    z = m.a*z + m.b(wh, qh)*m.tau2u(wh, tau(t(k)));
end
end
