% =========================================================================
% test_actuator_energy.m
% Metric 3 Evaluation: True mechanical actuator energy consumption
% (integral of sum |M_{m,i} * p_i| dt) across the 50% to 150% sweep.
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
scales   = [0.50, 0.75, 1.00, 1.25, 1.50];
n_runs   = length(scales);

total_energy_J  = zeros(n_runs, 1);
wheel_energy_J  = zeros(n_runs, 4);
effort_rw_N2m2s = zeros(n_runs, 1);

fprintf('=== Metric 3 Sweep: Actuator Mechanical Energy (50%% to 150%%) ===\n');

for i = 1:n_runs
    s = scales(i);
    fprintf('Testing scale %.2fx ... ', s);

    out_rw = run_sim_scenario(s, s, 'rw', mdl_name);
    [total_energy_J(i), wheel_energy_J(i, :), ~, effort_rw_N2m2s(i)] = ...
        calc_actuator_energy(out_rw.t, out_rw.M_m, out_rw.p, out_rw.u);

    fprintf('Done! [Total RW Energy: %.2f J]\n', total_energy_J(i));
end

T_energy = table(scales(:), total_energy_J, wheel_energy_J(:,1), wheel_energy_J(:,2), ...
    wheel_energy_J(:,3), wheel_energy_J(:,4), effort_rw_N2m2s, ...
    'VariableNames', {'Scale', 'Total_Energy_J', 'RW1_J', 'RW2_J', 'RW3_J', 'RW4_J', 'Control_Effort_u2'});
fprintf('\n');
disp(T_energy);

if show_plots
    figure('Name', 'Metric 3: Actuator Energy Sweep');
    plot(scales*100, total_energy_J, 'rs-', 'LineWidth', 1.8); grid on;
    xlabel('Target Inertia & Disturbance Scale (%)');
    ylabel('Total Mechanical Energy (Joules)');
    title('Metric 3: Tetrahedral Reaction Wheel Energy Consumption');
end