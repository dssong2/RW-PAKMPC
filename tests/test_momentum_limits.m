% =========================================================================
% test_momentum_limits.m
% Metric 4 Evaluation: Reaction wheel angular momentum tracking and
% saturation verification (|h_{w,i}| <= h_max) across the 50% to 150% sweep.
% =========================================================================

clear; clc;
show_plots = true;

test_dir = fileparts(mfilename('fullpath'));
if isempty(test_dir), test_dir = pwd; end
root_dir = fileparts(test_dir);
addpath(root_dir);
addpath(fullfile(root_dir, 'utils'));

evalin('base', 'config;');
mdl_name = 'pakmpc_rw_combined';
hmax_val = evalin('base', 'hmax');
scales   = [0.50, 0.75, 1.00, 1.25, 1.50];
n_runs   = length(scales);

peak_ratio_pct = zeros(n_runs, 1);
max_hw_wheels  = zeros(n_runs, 4);
saturated_flag = false(n_runs, 1);
sat_duration_s = zeros(n_runs, 1);

fprintf('=== Metric 4 Sweep: Reaction Wheel Momentum Management (h_max = %.1f Nms) ===\n', hmax_val);

for i = 1:n_runs
    s = scales(i);
    fprintf('Testing scale %.2fx ... ', s);

    out_rw = run_sim_scenario(s, s, 'rw', mdl_name);
    [p_ratio, max_w, is_sat, t_sat] = check_momentum_limits(out_rw.t, out_rw.hw, hmax_val);

    peak_ratio_pct(i)  = p_ratio * 100;
    max_hw_wheels(i,:) = max_w;
    saturated_flag(i)  = is_sat;
    sat_duration_s(i)  = t_sat;

    fprintf('Done! [Peak |h_w|/h_max: %.2f%% | Saturated: %d]\n', peak_ratio_pct(i), is_sat);
end

T_momentum = table(scales(:), peak_ratio_pct, max_hw_wheels(:,1), max_hw_wheels(:,2), ...
    max_hw_wheels(:,3), max_hw_wheels(:,4), saturated_flag, sat_duration_s, ...
    'VariableNames', {'Scale', 'Peak_h_Ratio_Pct', 'Max_RW1_Nms', 'Max_RW2_Nms', ...
    'Max_RW3_Nms', 'Max_RW4_Nms', 'Is_Saturated', 'Sat_Time_s'});
fprintf('\n');
disp(T_momentum);

if show_plots
    figure('Name', 'Metric 4: Momentum Management Sweep');
    plot(scales*100, peak_ratio_pct, 'bo-', 'LineWidth', 1.8); hold on;
    yline(100, 'r--', 'Saturation Limit (100%)', 'LineWidth', 1.5); grid on;
    xlabel('Target Inertia & Disturbance Scale (%)');
    ylabel('Peak Wheel Momentum (% of h_{max})');
    title('Metric 4: Peak Reaction Wheel Momentum vs. Saturation Limit');
end