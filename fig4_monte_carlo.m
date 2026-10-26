%FIG4_MONTE_CARLO  Fig. 4 of Sun et al.: settling time and total control effort of non-adaptive
%   Koopman MPC with the RBF-EDMD, DCM-based and quaternion-based models over 50 Monte Carlo runs
%   from random attitudes and random |w_i| <= 0.02 rad/s, no uncertainty (Jc0 = Jc).
%   Writes results/fig4_monte_carlo.png and results/fig4_monte_carlo.mat.
%
%   [paper] = quoted from the paper
%   [guess] = not given in the paper, chosen here

here = fileparts(mfilename('fullpath'));
if isempty(here) || contains(here, 'Editor_'), here = pwd; end
addpath(fullfile(here, 'lib'));

p = pakmpc_params();
nrun = 50;                                      % [paper]
rs = RandStream('mt19937ar', 'Seed', 42);       % [guess] seed
x0 = zeros(7, nrun);
for i = 1:nrun
    q = randn(rs, 4, 1); q = q/norm(q);
    if q(1) < 0, q = -q; end                    % [guess] shortest-rotation sign
    x0(:, i) = [0.02*(2*rand(rs, 3, 1) - 1); q];
end

models = {rbf_edmd(p), dcm_koopman(5, p.jc, p.ts), quat_koopman(5, p.jc, p.ts)};
names = {'RBF-EDMD', 'DCM-based', 'Quaternion-based'};
tset = nan(nrun, numel(models)); effort = nan(nrun, numel(models));
for j = 1:numel(models)
    c = koopman_mpc(models{j}, p);
    parfor i = 1:nrun
        o = simulate_attitude(c, x0(:, i), p);
        tset(i, j) = settling_time(o.t, o.x(5:7, :));
        effort(i, j) = sum(vecnorm(o.tau).^2)*p.ts;    % [paper] int ||tau||^2 dt
    end
    fprintf('%-17s settling median %.1f s [%.1f %.1f], unsettled %d; effort median %.3g [%.3g %.3g]\n', ...
        names{j}, median(tset(:, j), 'omitnan'), min(tset(:, j)), max(tset(:, j)), sum(isnan(tset(:, j))), ...
        median(effort(:, j)), min(effort(:, j)), max(effort(:, j)));
end

out = fullfile(here, 'results');
if ~exist(out, 'dir'), mkdir(out); end
save(fullfile(out, 'fig4_monte_carlo.mat'), 'x0', 'tset', 'effort', 'names');

fig = figure('Position', [100 100 720 330], 'Visible', 'off');
theme(fig, 'light');
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact');
g = categorical(repelem(1:numel(names), nrun)', 1:numel(names), names);
ax = nexttile(tl);
boxchart(ax, g, tset(:), 'MarkerStyle', '+');
ylabel(ax, 'Settling times (s)');
ax = nexttile(tl);
boxchart(ax, g, effort(:), 'MarkerStyle', '+');
ylabel(ax, 'Total Control Effort (N^2\cdotm^2\cdots)');
exportgraphics(fig, fullfile(out, 'fig4_monte_carlo.png'), 'Resolution', 200);
fprintf('saved %s\n', fullfile(out, 'fig4_monte_carlo.png'));
