% =========================================================================
% test_nominal_sanity.m
% Runs a single 1.0x nominal evaluation comparing Idealized Direct Torque vs.
% Realistic Tetrahedral Reaction Wheels and prints all 4 evaluation metrics.
% =========================================================================

clear; clc;
show_plots = true; % Set to true for this sanity check to visually verify metrics

% Setup paths automatically
test_dir = fileparts(mfilename('fullpath'));
if isempty(test_dir), test_dir = pwd; end
root_dir = fileparts(test_dir);
addpath(root_dir);
addpath(fullfile(root_dir, 'utils'));

evalin('base', 'config;');
mdl_name   = 'pakmpc_rw_combined';
settle_tol = 0.05;
hmax_val   = evalin('base', 'hmax');

fprintf('=== Running Nominal 1.0x Sanity Check (Ideal vs. Reaction Wheels) ===\n');

% 1. Run Idealized Baseline & Realistic RW Model
out_ideal = run_sim_scenario(1.0, 1.0, 'ideal', mdl_name);
out_rw    = run_sim_scenario(1.0, 1.0, 'rw',    mdl_name);

% 2. Compute Metrics via utils/
[ss_id, tconv_id, qn_id] = calc_convergence(out_ideal.t, out_ideal.q, settle_tol);
[ss_rw, tconv_rw, qn_rw] = calc_convergence(out_rw.t,    out_rw.q,    settle_tol);

[energy_rw, ~, ~, effort_rw] = calc_actuator_energy(out_rw.t, out_rw.M_m, out_rw.p, out_rw.u);
[~, ~, ~, effort_id]         = calc_actuator_energy(out_ideal.t, out_ideal.M_m, out_ideal.p, out_ideal.u);

[peak_h_ratio, ~, is_sat, sat_t] = check_momentum_limits(out_rw.t, out_rw.hw, hmax_val);

% 3. Display Summary Table
fprintf('\n--- Nominal Sanity Check Summary ---\n');
fprintf('Metric                               | Ideal Baseline     | Realistic RWs\n');
fprintf('--------------------------------------------------------------------------\n');
fprintf('1. Steady-State Error ||q_v||_2      | %.3e          | %.3e\n', ss_id, ss_rw);
fprintf('1. Convergence Time (5%% threshold)   | %.2f s             | %.2f s\n', tconv_id, tconv_rw);
fprintf('2. Avg Solver Step Runtime           | %.4f s           | %.4f s\n', out_ideal.avg_step_time, out_rw.avg_step_time);
fprintf('3. Classical Effort (int ||u||^2 dt) | %.3e N^2m^2s    | %.3e N^2m^2s\n', effort_id, effort_rw);
fprintf('3. Mechanical RW Energy              | N/A (Black-Box)    | %.3e J\n', energy_rw);
fprintf('4. Peak Wheel Momentum |h_w|/h_max   | N/A                | %.2f%% (Sat: %d, %.1fs)\n', peak_h_ratio*100, is_sat, sat_t);

% 4. Optional Visual Verification
if show_plots
    figure('Name', 'Convergence Sanity Check Overlay', 'Position', [100, 100, 750, 420]);
    plot(out_ideal.t, qn_id, 'b--', 'LineWidth', 1.6); hold on;
    plot(out_rw.t,    qn_rw, 'r-',  'LineWidth', 1.6); grid on;
    yline(settle_tol * qn_rw(1), 'k:', '5% Threshold', 'LineWidth', 1.4);
    if ~isnan(tconv_id), xline(tconv_id, 'b:', sprintf('t_{id}=%.1fs', tconv_id), 'LineWidth', 1.2); end
    if ~isnan(tconv_rw), xline(tconv_rw, 'r:', sprintf('t_{rw}=%.1fs', tconv_rw), 'LineWidth', 1.2); end
    xlabel('Time (s)'); ylabel('Attitude Error ||q_v||_2');
    title('Sanity Check: Attitude Error & Convergence Thresholds');
    legend('Ideal Direct Torque', 'Tetrahedral RWs', 'Location', 'northeast');

    plot_telemetry(out_rw, 'Nominal 1.0x RW');
end