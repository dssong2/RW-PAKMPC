function plot_telemetry(sim_out, title_prefix)
% PLOT_TELEMETRY Opt-in visualization helper for a single simulation run or pair.
%
% Usage:
%   plot_telemetry(out_rw, 'Nominal RW Run');

if nargin < 2 || isempty(title_prefix)
    title_prefix = sprintf('Telemetry (%s mode, Scale=%.2fx)', upper(sim_out.mode), sim_out.target_scale);
end

taumax = evalin('base', 'taumax');
hmax   = evalin('base', 'hmax');

figure('Name', title_prefix, 'Position', [120, 100, 850, 720]);

subplot(4, 1, 1);
plot(sim_out.t, sim_out.q(:, 2:4), 'LineWidth', 1.4); grid on;
ylabel('q_v');
title([title_prefix ' - Attitude & Reaction Wheel Telemetry']);
legend('q_1', 'q_2', 'q_3', 'Location', 'northeast');

subplot(4, 1, 2);
plot(sim_out.t, sim_out.w * 180/pi, 'LineWidth', 1.4); grid on;
ylabel('Body \omega (deg/s)');
legend('\omega_x', '\omega_y', '\omega_z', 'Location', 'northeast');

subplot(4, 1, 3);
plot(sim_out.t, sim_out.hw, 'LineWidth', 1.4); hold on;
yline([hmax, -hmax], 'k--', 'h_{max}'); grid on;
ylabel('h_w (N\cdotm\cdots)');
legend('RW 1', 'RW 2', 'RW 3', 'RW 4', 'Location', 'northeast');

subplot(4, 1, 4);
stairs(sim_out.tm, sim_out.M_m_d, 'LineWidth', 1.4); hold on;
yline([taumax, -taumax], 'k--', '\tau_{max}'); grid on;
ylabel('Motor M_m (N\cdotm)'); xlabel('Time (s)');
legend('RW 1', 'RW 2', 'RW 3', 'RW 4', 'Location', 'northeast');
end