%FIG6_11_TABLE3_ADAPTIVE  Figs. 6-11 and Table III of Sun et al.: PAKMPC against non-adaptive MPC
%   on the combined spacecraft with inertia uncertainty (Jc0 = Js, plant Jc), the target torque of
%   Eq. (59), the disturbance tau_d and measurement noise, Sec. V-C.
%   Writes results/fig6_qv.png ... results/fig11_db.png and prints Table III.
%   KMAE-AKMPC and INMPC (step 4) are added to the controllers list when they exist.
%
%   [paper] = quoted from the paper
%   [guess] = not given in the paper, chosen here

here = fileparts(mfilename('fullpath'));
if isempty(here) || contains(here, 'Editor_'), here = pwd; end
addpath(fullfile(here, 'lib'));

p = pakmpc_params();
x0 = [p.w0; p.q0];
m = quat_koopman(p.nphi, p.jc0, p.ts);          % [paper] Nphi = 5, Jc0 = Js
plant = struct('dist', p.dist, 'target', true, 'noise', true, ...
    'target_ref', 'identity');                  % [guess] target regulates to the identity: this
                                                %   matches Table III's steady-state ||q_v||; 'initial'
                                                %   (build_pakmpc_multibody) leaves ~7e-3
ssw = [190 200];                                % [guess] steady-state window, as in the Fig. 6 insets

ctrl = {'Proposed PAKMPC', koopman_mpc(m, p, struct('adapt', true)), '-', [0.85 0.2 0.1];
        'Non-adaptive MPC', koopman_mpc(m, p), '--', [0.2 0.4 0.85]};
nc = size(ctrl, 1);
runs = cell(1, nc); rt = zeros(1, nc);
for i = 1:nc
    t0 = tic;
    runs{i} = simulate_attitude(ctrl{i, 2}, x0, p, plant);
    rt(i) = toc(t0);
end
out = fullfile(here, 'results');
if ~exist(out, 'dir'), mkdir(out); end

% ---------------- Table III ----------------
fprintf('\nTABLE III\n%-18s %12s %10s %12s %12s %12s\n', 'Strategy', 'runtime (s)', 'step (s)', ...
    'ss ||q_v||', 'ss ||w||', 'effort');
for i = 1:nc
    o = runs{i};
    s = o.t >= ssw(1) & o.t <= ssw(2);
    fprintf('%-18s %12.2f %10.4f %12.2e %12.2e %12.2e\n', ctrl{i, 1}, rt(i), mean(o.steptime), ...
        mean(vecnorm(o.x(5:7, s))), mean(vecnorm(o.x(1:3, s))), sum(vecnorm(o.tau).^2)*p.ts);
end

% ---------------- Figs. 6-8: q_v, w, tau with steady-state insets ----------------
specs = {'fig6_qv', @(o) o.x(5:7, :), @(o) o.t, {'q_1', 'q_2', 'q_3'}, true;
         'fig7_w', @(o) o.x(1:3, :), @(o) o.t, {'\omega_1 (rad/s)', '\omega_2 (rad/s)', '\omega_3 (rad/s)'}, true;
         'fig8_tau', @(o) o.tau, @(o) o.t(1:end-1), {'\tau_1 (N\cdotm)', '\tau_2 (N\cdotm)', '\tau_3 (N\cdotm)'}, false};
for f = 1:size(specs, 1)
    fig = figure('Position', [100 100 560 620], 'Visible', 'off');
    theme(fig, 'light');
    tl = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact');
    axs = gobjects(1, 3);
    for j = 1:3
        ax = nexttile(tl); hold(ax, 'on'); box(ax, 'on'); axs(j) = ax;
        for i = 1:nc
            y = specs{f, 2}(runs{i}); t = specs{f, 3}(runs{i});
            plot(ax, t, y(j, :), ctrl{i, 3}, 'Color', ctrl{i, 4}, 'LineWidth', 1.2, 'DisplayName', ctrl{i, 1});
        end
        ylabel(ax, specs{f, 4}{j}); xlim(ax, [0 p.tend]);
        if strcmp(specs{f, 1}, 'fig8_tau'), ylim(ax, [-20 20]); end
    end
    xlabel(axs(3), 'Time (s)');
    lg = legend(axs(1), 'NumColumns', 2); lg.Layout.Tile = 'north';
    if specs{f, 5}
        drawnow;
        for j = 1:3
            pos = axs(j).Position;
            ins = axes(fig, 'Position', [pos(1) + 0.5*pos(3), pos(2) + 0.45*pos(4), 0.42*pos(3), 0.45*pos(4)]);
            hold(ins, 'on'); box(ins, 'on');
            for i = 1:nc
                y = specs{f, 2}(runs{i}); t = specs{f, 3}(runs{i});
                s = t >= ssw(1);
                plot(ins, t(s), y(j, s), ctrl{i, 3}, 'Color', ctrl{i, 4}, 'LineWidth', 1);
            end
            xlim(ins, ssw); ins.FontSize = 7;
        end
    end
    exportgraphics(fig, fullfile(out, [specs{f, 1} '.png']), 'Resolution', 200);
end

% ---------------- Fig. 9: body-axis trajectories on the unit sphere ----------------
fig = figure('Position', [100 100 760 300], 'Visible', 'off');
theme(fig, 'light');
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact');
[sx, sy, sz] = sphere(24);
axn = {'X-axis', 'Y-axis', 'Z-axis'};
for a = 1:3
    ax = nexttile(tl); hold(ax, 'on');
    surf(ax, sx, sy, sz, 'FaceAlpha', 0.05, 'EdgeColor', [0.8 0.8 0.8], 'FaceColor', [0.9 0.9 0.9], ...
        'HandleVisibility', 'off');
    for i = 1:nc
        o = runs{i};
        v = zeros(3, numel(o.t));
        for k = 1:numel(o.t)
            r = q2r(o.x(4:7, k)); v(:, k) = r(:, a);
        end
        plot3(ax, v(1, :), v(2, :), v(3, :), ctrl{i, 3}, 'Color', ctrl{i, 4}, 'LineWidth', 1.5, ...
            'DisplayName', ctrl{i, 1});
    end
    r0 = q2r(x0(4:7));
    plot3(ax, r0(1, a), r0(2, a), r0(3, a), 'p', 'Color', [0.85 0.3 0.1], 'MarkerSize', 10, ...
        'MarkerFaceColor', [0.85 0.3 0.1], 'DisplayName', 'Initial Attitude');
    axis(ax, 'equal', 'off'); view(ax, 135, 25); title(ax, axn{a}, 'FontWeight', 'normal');
end
lg = legend(nexttile(tl, 1), 'NumColumns', 3); lg.Layout.Tile = 'north';
exportgraphics(fig, fullfile(out, 'fig9_sphere.png'), 'Resolution', 200);

% ---------------- Figs. 10-11: mu_k and ||Delta B|| of PAKMPC ----------------
o = runs{1};
names = {'fig10_mu', o.mu, '\mu_k'; 'fig11_db', o.db, '||\DeltaB||_2'};
for f = 1:2
    fig = figure('Position', [100 100 520 220], 'Visible', 'off');
    theme(fig, 'light');
    ax = axes(fig); box(ax, 'on');
    plot(ax, o.t(1:end-1), names{f, 2}, '-', 'Color', ctrl{1, 4}, 'LineWidth', 1.2);
    xlabel(ax, 'Time (s)'); ylabel(ax, names{f, 3}); xlim(ax, [0 p.tend]);
    exportgraphics(fig, fullfile(out, [names{f, 1} '.png']), 'Resolution', 200);
end
fprintf('\nsaved fig6_qv, fig7_w, fig8_tau, fig9_sphere, fig10_mu, fig11_db in %s\n', out);
