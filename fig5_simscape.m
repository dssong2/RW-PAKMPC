%FIG5_SIMSCAPE  PAKMPC on the Simscape Multibody plant of build_pakmpc_multibody, Sec. V-C and Fig. 5.
%   Loads pakmpc_multibody.slx (built by build_pakmpc_multibody if missing), puts PAKMPC in its
%   Controller block through simscape_step, matches the plant settings to fig6_11_table3_adaptive,
%   and runs the closed loop. Checks it against the ode45 plant of simulate_attitude without noise
%   (same inputs, so the two should agree to integration accuracy), then runs the noisy case and
%   draws the Fig. 5 snapshots at t = 0, 10, 20, 50 s.
%   Writes results/fig5_snapshots.png and results/simscape_vs_ode45.png. The model is changed in
%   memory only; pakmpc_multibody.slx on disk keeps Elijah's PD controller.
%
%   [paper] = quoted from the paper
%   [guess] = not given in the paper, chosen here

here = fileparts(mfilename('fullpath'));
if isempty(here) || contains(here, 'Editor_'), here = pwd; end
addpath(fullfile(here, 'lib'));

p = pakmpc_params();
x0 = [p.w0; p.q0];
m = quat_koopman(p.nphi, p.jc0, p.ts);
mdl = 'pakmpc_multibody';
if ~exist(fullfile(here, [mdl '.slx']), 'file')
    run(fullfile(here, 'build_pakmpc_multibody.m'));
end
load_system(fullfile(here, [mdl '.slx']));
set_param(mdl, 'SimMechanicsOpenEditorOnUpdate', 'off');
setscript([mdl '/Controller'], {
    'function tau = ctl(wm, qm)'
    '% PAKMPC through lib/simscape_step'
    'coder.extrinsic(''simscape_step'');'
    'tau = zeros(3, 1);'
    'tau = simscape_step(wm, qm);'});
% [guess] same plant settings as fig6_11_table3_adaptive: target regulates to the identity and
% noise std 1e-4 (pakmpc_params); build_pakmpc_multibody holds the initial attitude, variance 1e-4
setscript([mdl '/Target Torque'], {
    'function ut = target(q, w)'
    '% Eq. (59), u_t = -Kpt q_vt - Kdt w, target regulates to the identity'
    'q = q(:); w = w(:);'
    ['kpt = ' mat2str(p.kpt, 17) '; kdt = ' mat2str(p.kdt, 17) ';']
    'ut = -kpt*q(2:4) - kdt*w;'});
set_param([mdl '/Noise'], 'Variance', mat2str(p.noise_var*ones(1, 6)));

% ---------------- check: Simscape against ode45, no noise ----------------
set_param([mdl '/Noise On'], 'Gain', '0');
simscape_step('init', koopman_mpc(m, p, struct('adapt', true)));
so = sim(mdl, 'StopTime', num2str(p.tend));
ss0 = unpack(so, simscape_step('log'));
oo = simulate_attitude(koopman_mpc(m, p, struct('adapt', true)), x0, p, ...
    struct('dist', p.dist, 'target', true, 'target_ref', 'identity'));
dq = max(vecnorm(ss0.x(4:7, :) - oo.x(4:7, :)));
dtau = max(vecnorm(ss0.tau - oo.tau));
fprintf('no noise: max |q_simscape - q_ode45| = %.2e, max |tau diff| = %.2e N m\n', dq, dtau);
fprintf('          ss ||q_v||: Simscape %.2e, ode45 %.2e\n', ssq(ss0), ssq(oo));

% ---------------- noisy run (Simscape noise stream differs from simulate_attitude) ----------------
set_param([mdl '/Noise On'], 'Gain', '1');
simscape_step('init', koopman_mpc(m, p, struct('adapt', true)));
so = sim(mdl, 'StopTime', num2str(p.tend));
ss = unpack(so, simscape_step('log'));
on = simulate_attitude(koopman_mpc(m, p, struct('adapt', true)), x0, p, ...
    struct('dist', p.dist, 'target', true, 'noise', true, 'target_ref', 'identity'));
fprintf('noise:    ss ||q_v||: Simscape %.2e, ode45 %.2e; ss ||w||: Simscape %.2e, ode45 %.2e\n', ...
    ssq(ss), ssq(on), ssw(ss), ssw(on));
fprintf('          mean step %.4f s, effort Simscape %.3g, ode45 %.3g\n', mean(ss.steptime), ...
    sum(vecnorm(ss.tau).^2)*p.ts, sum(vecnorm(on.tau).^2)*p.ts);
close_system(mdl, 0);

out = fullfile(here, 'results');
if ~exist(out, 'dir'), mkdir(out); end

% ---------------- Simscape against ode45 ----------------
fig = figure('Position', [100 100 560 520], 'Visible', 'off');
theme(fig, 'light');
tl = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact');
sets = {ss.x(5:7, :), on.x(5:7, :), 'q_v'; ss.x(1:3, :), on.x(1:3, :), '\omega (rad/s)'; ss.tau, on.tau, '\tau (N\cdotm)'};
for j = 1:3
    ax = nexttile(tl); hold(ax, 'on'); box(ax, 'on');
    t = ss.t(1:size(sets{j, 1}, 2));
    plot(ax, t, sets{j, 1}, '-', 'LineWidth', 1.2);
    ax.ColorOrderIndex = 1;
    plot(ax, t, sets{j, 2}, '--', 'LineWidth', 1.2);
    ylabel(ax, sets{j, 3}); xlim(ax, [0 80]);
end
xlabel(tl, 'Time (s)'); title(tl, 'PAKMPC with noise: Simscape (solid), ode45 (dashed)');
exportgraphics(fig, fullfile(out, 'simscape_vs_ode45.png'), 'Resolution', 200);

% ---------------- Fig. 5: snapshots ----------------
fig = figure('Position', [100 100 1000 290], 'Visible', 'off', 'Color', 'w');
theme(fig, 'light');
tl = tiledlayout(fig, 1, 4, 'TileSpacing', 'none', 'Padding', 'compact');
ts5 = [0 10 20 50];                             % [paper] Fig. 5
lab = 'abcd';
for i = 1:4
    ax = nexttile(tl);
    [~, k] = min(abs(ss.t - ts5(i)));
    draw_combined(ax, ss.x(4:7, k), p);
    axis(ax, 'equal', 'off'); axis(ax, 2.9*[-1 1 -1 1 -1 1]);
    view(ax, 135, 25); camlight(ax, 'headlight'); lighting(ax, 'gouraud');
    title(ax, sprintf('(%s) t = %d s', lab(i), ts5(i)), 'FontWeight', 'normal');
end
exportgraphics(fig, fullfile(out, 'fig5_snapshots.png'), 'Resolution', 200);
fprintf('\nsaved fig5_snapshots, simscape_vs_ode45 in %s\n', out);

function o = unpack(so, lg)
% sample the logged Simscape signals at the controller instants, as simulate_attitude returns them
t = 0:0.1:so.q_log.time(end);
q = interp1(so.q_log.time, sig(so.q_log, 4), t)';
w = interp1(so.w_log.time, sig(so.w_log, 3), t)';
if size(lg.tau, 2) ~= numel(t) - 1 && size(lg.tau, 2) ~= numel(t)
    warning('controller ran %d times for %d samples', size(lg.tau, 2), numel(t) - 1);
end
o = struct('t', t, 'x', [w; q], 'tau', lg.tau(:, 1:numel(t) - 1), 'mu', lg.mu, 'db', lg.db, 'steptime', lg.steptime);
end

function v = ssq(o), s = o.t >= 190; v = mean(vecnorm(o.x(5:7, s))); end
function v = ssw(o), s = o.t >= 190; v = mean(vecnorm(o.x(1:3, s))); end

function x = sig(s, n)
x = squeeze(s.signals.values);
if size(x, 2) ~= n, x = x.'; end
end

function setscript(path, lines)
ch = sfroot().find('-isa', 'Stateflow.EMChart', 'Path', path);
ch.Script = strjoin(lines, newline);
end
