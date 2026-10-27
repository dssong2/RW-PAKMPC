%FIG12_13_ROBUSTNESS  Figs. 12-13 of Sun et al., Sec. V-D: structured sweep of the target inertia
%   and mass (jointly scaled 50-150%) and the disturbance amplitude (0-200%) for PAKMPC, non-adaptive
%   MPC, KMAE-AKMPC and INMPC, with everything else as in fig6_11_table3_adaptive.
%   Writes results/fig12a_heatmap_qv.png, results/fig12b_heatmap_effort.png,
%   results/fig13_scatter.png and results/fig12_13_robustness.mat.
%
%   [paper] = quoted from the paper
%   [guess] = not given in the paper, chosen here

here = fileparts(mfilename('fullpath'));
if isempty(here) || contains(here, 'Editor_'), here = pwd; end
addpath(fullfile(here, 'lib'));

p = pakmpc_params();
x0 = [p.w0; p.q0];
m = quat_koopman(p.nphi, p.jc0, p.ts);
si = [0.5 0.75 1 1.25 1.5];                     % [paper] target inertia and mass scaling
sd = [0 0.5 1 1.5 2];                           % [paper] disturbance scaling
ssw = [190 200];                                % [guess] steady-state window, as in fig6_11
names = {'PAKMPC', 'Non-adaptive MPC', 'KMAE-AKMPC', 'INMPC'};
ctrls = {koopman_mpc(m, p, struct('adapt', true)), koopman_mpc(m, p), kmae_akmpc(p), incremental_nmpc(p)};
colors = [0.85 0.33 0.1; 0.2 0.4 0.85; 0.55 0.25 0.65; 0.4 0.75 0.3];

[ii, jj, kk] = ndgrid(1:numel(si), 1:numel(sd), 1:numel(ctrls));
nr = numel(ii);
qv = zeros(1, nr); eff = zeros(1, nr); st = zeros(1, nr);
parfor r = 1:nr
    ps = p;
    [jc, jt] = combined_inertia(p, si(ii(r)));
    ps.kpt = 0.02*jt; ps.kdt = 0.05*jt;          % [paper] Eq. (59) gains follow the scaled Jt
    a = sd(jj(r));
    plant = struct('jc', jc, 'dist', @(t) a*p.dist(t), 'target', true, 'noise', true, ...
        'target_ref', 'identity');
    o = simulate_attitude(ctrls{kk(r)}, x0, ps, plant);
    s = o.t >= ssw(1) & o.t <= ssw(2);
    qv(r) = mean(vecnorm(o.x(5:7, s)));
    eff(r) = sum(vecnorm(o.tau).^2)*p.ts;
    st(r) = mean(o.steptime);
end
qv = reshape(qv, size(ii)); eff = reshape(eff, size(ii)); st = reshape(st, size(ii));
out = fullfile(here, 'results');
if ~exist(out, 'dir'), mkdir(out); end
save(fullfile(out, 'fig12_13_robustness.mat'), 'si', 'sd', 'names', 'qv', 'eff', 'st');

fprintf('\n%-18s %12s %12s %12s %12s\n', 'Strategy', 'mean ||q_v||', 'max ||q_v||', 'mean effort', 'step (ms)');
for k = 1:numel(names)
    a = qv(:, :, k); b = eff(:, :, k); c = st(:, :, k);
    fprintf('%-18s %12.2e %12.2e %12.2e %12.2f\n', names{k}, mean(a(:)), max(a(:)), mean(b(:)), 1e3*mean(c(:)));
end

% ---------------- Fig. 12: heatmaps, shared color scale per metric ----------------
mets = {'fig12a_heatmap_qv', qv, 'Steady-state error ||q_v||';
        'fig12b_heatmap_effort', eff, 'Total control effort \int||\tau||_2^2 dt'};
for f = 1:2
    v = mets{f, 2};
    fig = figure('Position', [100 100 640 330], 'Visible', 'off');
    theme(fig, 'light');
    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact');
    for k = 1:numel(names)
        ax = nexttile(tl);
        imagesc(ax, 100*si, 100*sd, v(:, :, k)');
        set(ax, 'YDir', 'normal', 'XTick', 100*si, 'YTick', 100*sd);
        clim(ax, [min(v(:)) max(v(:))]); colormap(ax, parula);
        title(ax, names{k}, 'FontWeight', 'normal');
    end
    cb = colorbar(ax); cb.Layout.Tile = 'east'; cb.Label.String = mets{f, 3};
    xlabel(tl, 'Target inertia scaling factor (%)'); ylabel(tl, 'Disturbance scaling factor (%)');
    exportgraphics(fig, fullfile(out, [mets{f, 1} '.png']), 'Resolution', 200);
end

% ---------------- Fig. 13: runtime against effort and accuracy ----------------
fig = figure('Position', [100 100 640 300], 'Visible', 'off');
theme(fig, 'light');
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact');
ys = {eff, 'Total control effort (N^2m^2s)', 'linear', 'Runtime vs. control effort';
      qv, 'Steady-state error ||q_v||', 'log', 'Runtime vs. attitude precision'};
for f = 1:2
    ax = nexttile(tl); hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
    for k = 1:numel(names)
        x = 1e3*st(:, :, k); y = ys{f, 1}(:, :, k);
        scatter(ax, x(:), y(:), 15, colors(k, :), 'filled', 'MarkerFaceAlpha', 0.3, 'HandleVisibility', 'off');
        scatter(ax, mean(x(:)), mean(y(:)), 110, colors(k, :), 'filled', 'MarkerEdgeColor', 'k', ...
            'LineWidth', 1.2, 'DisplayName', names{k});
    end
    set(ax, 'YScale', ys{f, 3});
    xlabel(ax, 'Average step runtime (ms)'); ylabel(ax, ys{f, 2}); title(ax, ys{f, 4}, 'FontWeight', 'normal');
end
lg = legend(ax, 'NumColumns', 4); lg.Layout.Tile = 'north';
exportgraphics(fig, fullfile(out, 'fig13_scatter.png'), 'Resolution', 200);
fprintf('\nsaved fig12a_heatmap_qv, fig12b_heatmap_effort, fig13_scatter in %s\n', out);
