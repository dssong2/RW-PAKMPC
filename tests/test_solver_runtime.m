% =========================================================================
% test_solver_runtime.m
% Metric 2 Evaluation: Average step runtime across the 50% to 150% sweep.
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

step_time_ideal = zeros(n_runs, 1);
step_time_rw    = zeros(n_runs, 1);

% Warm-up run so JIT compilation time is not counted in the benchmark
fprintf('Warming up Simulink model...\n');
run_sim_scenario(1.0, 1.0, 'rw', mdl_name);

fprintf('=== Metric 2 Sweep: Solver Step Runtime (50%% to 150%%) ===\n');

for i = 1:n_runs
    s = scales(i);
    fprintf('Timing scale %.2fx ... ', s);

    out_id = run_sim_scenario(s, s, 'ideal', mdl_name);
    out_rw = run_sim_scenario(s, s, 'rw',    mdl_name);

    step_time_ideal(i) = out_id.avg_step_time * 1000; % Convert to ms
    step_time_rw(i)    = out_rw.avg_step_time * 1000; % Convert to ms

    fprintf('Done! [Ideal: %.3f ms/step | RW: %.3f ms/step]\n', ...
        step_time_ideal(i), step_time_rw(i));
end

T_runtime = table(scales(:), step_time_ideal, step_time_rw, ...
    'VariableNames', {'Scale', 'Ideal_Avg_Step_ms', 'RW_Avg_Step_ms'});
fprintf('\n');
disp(T_runtime);
fprintf('Overall Mean Step Time -> Ideal: %.3f ms | RW: %.3f ms\n', ...
    mean(step_time_ideal), mean(step_time_rw));

if show_plots
    figure('Name', 'Metric 2: Step Runtime Sweep');
    bar(scales*100, [step_time_ideal, step_time_rw]); grid on;
    xlabel('Target Inertia & Disturbance Scale (%)');
    ylabel('Average Step Runtime (ms)');
    title('Metric 2: Computational Runtime per Discrete Step (T_s = 100 ms)');
    legend('Ideal Direct Torque', 'Tetrahedral RWs', 'Location', 'northwest');
end