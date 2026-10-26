%FIG3_TABLE2_NOMINAL_MPC  Fig. 3 and Table II of Sun et al.: non-adaptive Koopman MPC with the
%   RBF-EDMD, DCM-based and quaternion-based nominal models, no uncertainty (Jc0 = Jc), no
%   disturbance or noise, from q0 = [0.5 0.5 0.5 0.5], w0 = [0.02 -0.01 -0.01] rad/s.
%   Writes results/fig3_nominal_mpc.png and prints Table II (runtime, step time, peak memory).
%
%   [paper] = quoted from the paper
%   [guess] = not given in the paper, chosen here

here = fileparts(mfilename('fullpath'));
if isempty(here) || contains(here, 'Editor_'), here = pwd; end
addpath(fullfile(here, 'lib'));

p = pakmpc_params();
x0 = [p.w0; p.q0];
tplot = 150;                                    % [paper] Fig. 3 time axis
edmd = rbf_edmd(p);                             % [guess] data set, see rbf_edmd

% ---------------- Fig. 3: Nphi = 5 ----------------
models = {edmd, dcm_koopman(5, p.jc, p.ts), quat_koopman(5, p.jc, p.ts)};
runs = cell(1, numel(models));
for i = 1:numel(models)
    runs{i} = simulate_attitude(koopman_mpc(models{i}, p), x0, p);
end

titles = {sprintf('RBF-based EDMD\nKoopman model'), sprintf('DCM-based\nKoopman model'), ...
    sprintf('Quaternion-based\nKoopman model')};
styles = {'-', '--', ':'};
cols = [0.85 0.2 0.1; 0.2 0.4 0.85; 0.1 0.7 0.75];
ylab = {'q_v', '\omega (rad/s)', '\tau (N\cdotm)'};
ylims = {[-0.2 0.5], [-0.15 0.1], [-20 20]};       % [paper] limits, w widened to show the dip
fig = figure('Position', [100 100 760 560], 'Visible', 'off');
theme(fig, 'light');
tl = tiledlayout(fig, 3, 3, 'TileSpacing', 'compact');
for r = 1:3
    for i = 1:numel(models)
        ax = nexttile(tl, (r - 1)*3 + i); hold(ax, 'on'); box(ax, 'on');
        o = runs{i};
        switch r
            case 1, y = o.x(5:7, :); t = o.t;
            case 2, y = o.x(1:3, :); t = o.t;
            case 3, y = o.tau; t = o.t(1:end-1);
        end
        for j = 1:3
            plot(ax, t, y(j, :), styles{j}, 'Color', cols(j, :), 'LineWidth', 1.1, ...
                'DisplayName', sprintf('i = %d', j));
        end
        xlim(ax, [0 tplot]); ylim(ax, ylims{r}); ylabel(ax, ylab{r});
        if r == 1, title(ax, titles{i}, 'FontWeight', 'normal'); end
        if r == 3 && i == 2, xlabel(ax, 'Time (s)'); end
    end
end
lg = legend(nexttile(tl, 1), 'Orientation', 'horizontal');
lg.Layout.Tile = 'north';
out = fullfile(here, 'results');
if ~exist(out, 'dir'), mkdir(out); end
exportgraphics(fig, fullfile(out, 'fig3_nominal_mpc.png'), 'Resolution', 200);

% ---------------- Table II ----------------
% [paper] runtime of the full 200 s closed loop, mean controller step, peak memory (Profiler)
rows = {'Quaternion-based model', @(nphi) quat_koopman(nphi, p.jc, p.ts), [5 7];
        'DCM-based model [18]', @(nphi) dcm_koopman(nphi, p.jc, p.ts), [5 7];
        'RBF-based EDMD model [22]', @(nphi) edmd, nan};
fprintf('\nTABLE II\n%-27s %6s %8s %12s %10s %12s\n', 'Model', 'Nphi', 'lifted', 'runtime (s)', 'step (s)', 'memory (KB)');
for i = 1:size(rows, 1)
    for nphi = rows{i, 3}
        m = rows{i, 2}(nphi);
        c = koopman_mpc(m, p);
        simulate_attitude(c, x0, p, struct('tend', 1));          % warm up
        t0 = tic; o = simulate_attitude(c, x0, p); rt = toc(t0);
        profile('clear'); profile('-memory', 'on');     % separate run, the profiler slows it down
        simulate_attitude(c, x0, p, struct('tend', 20));
        profile('off');
        mem = peak_memory(profile('info'));
        fprintf('%-27s %6s %8d %12.2f %10.4f %12.0f\n', rows{i, 1}, num2str(nphi), m.n, rt, ...
            mean(o.steptime), mem);
    end
end
fprintf('\nsaved %s\n', fullfile(out, 'fig3_nominal_mpc.png'));

function kb = peak_memory(info)
% largest PeakMem over the profiled functions, KB; NaN if the profiler gave none
kb = nan;
ft = info.FunctionTable;
if ~isempty(ft) && isfield(ft, 'PeakMem')
    kb = max([ft.PeakMem])/1024;
end
end
