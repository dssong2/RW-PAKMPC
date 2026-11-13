% =========================================================================
% test_convergence.m
% Metric 1 Evaluation: Steady-state quaternion vector error and convergence
% time across a 50% to 150% parameter sweep for Ideal vs. RW actuation.
% =========================================================================

clear; clc;
show_plots = true; % Set true to plot sweep comparison curves

test_dir = fileparts(mfilename('fullpath'));
if isempty(test_dir), test_dir = pwd; end
root_dir = fileparts(test_dir);
addpath(root_dir);
addpath(fullfile(root_dir, 'utils'));

evalin('base', 'config;');
mdl_name   = 'pakmpc_rw_combined';
scales     = [0.50, 0.75, 1.00, 1.25, 1.50];
settle_tol = 0.05;

n_runs = length(scales);
ss_err_ideal = zeros(n_runs, 1); t_conv_ideal = zeros(n_runs, 1);
ss_err_rw    = zeros(n_runs, 1); t_conv_rw    = zeros(n_runs, 1);

fprintf('=== Metric 1 Sweep: Attitude Error & Convergence Time (50%% to 150%%) ===\n');

for i = 1:n_runs
    s = scales(i);
    fprintf('Running scale %.2fx ... ', s);

    out_id = run_sim_scenario(s, s, 'ideal', mdl_name);
    out_rw = run_sim_scenario(s, s, 'rw',    mdl_name);

    [ss_err_ideal(i), t_conv_ideal(i)] = calc_convergence(out_id.t, out_id.q, settle_tol);
    [ss_err_rw(i),    t_conv_rw(i)]    = calc_convergence(out_rw.t, out_rw.q, settle_tol);

    fprintf('Done! [Ideal: SS=%.2e, T=%.1fs | RW: SS=%.2e, T=%.1fs]\n', ...
        ss_err_ideal(i), t_conv_ideal(i), ss_err_rw(i), t_conv_rw(i));
end

T_results = table(scales(:), ss_err_ideal, t_conv_ideal, ss_err_rw, t_conv_rw, ...
    'VariableNames', {'Scale', 'Ideal_SS_Error', 'Ideal_T_Conv_s', 'RW_SS_Error', 'RW_T_Conv_s'});
fprintf('\n');
disp(T_results);

if show_plots
    figure('Name', 'Metric 1: Convergence Sweep', 'Position', [150, 150, 750, 500]);
    subplot(2,1,1);
    plot(scales*100, ss_err_ideal, 'bo--', 'LineWidth', 1.5); hold on;
    plot(scales*100, ss_err_rw,    'rs-',  'LineWidth', 1.5); grid on;
    ylabel('Steady-State ||q_v||_2'); title('Metric 1: Attitude Stabilization Across Parameter Sweep');
    legend('Ideal Direct Torque', 'Tetrahedral RWs', 'Location', 'northwest');

    subplot(2,1,2);
    plot(scales*100, t_conv_ideal, 'bo--', 'LineWidth', 1.5); hold on;
    plot(scales*100, t_conv_rw,    'rs-',  'LineWidth', 1.5); grid on;
    xlabel('Target Inertia & Disturbance Scale (%)'); ylabel('Settling Time (s)');
    legend('Ideal Direct Torque', 'Tetrahedral RWs', 'Location', 'northwest');
end