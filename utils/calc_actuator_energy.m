function [total_energy, wheel_energy, power_t, ideal_effort] = calc_actuator_energy(t, M_m, p, u)
% CALC_ACTUATOR_ENERGY Computes mechanical RW energy and idealized control effort.
%
% Inputs:
%   t   - (N x 1) Continuous simulation time vector (s)
%   M_m - (N x 4) Applied motor torques aligned to t (N*m)
%   p   - (N x 4) Wheel angular velocities (rad/s)
%   u   - (N x 3) [Optional] Commanded 3-axis body torque (N*m) for paper's effort metric
%
% Outputs:
%   total_energy - Total mechanical energy across all 4 wheels: \int \sum |M_m,i * p_i| dt (J)
%   wheel_energy - (1 x 4) Mechanical energy per wheel (J)
%   power_t      - (N x 1) Instantaneous mechanical power across all 4 wheels (W)
%   ideal_effort - Classical control effort \int ||u||_2^2 dt (N^2*m^2*s)

t = t(:);
num_steps = length(t);

if size(M_m, 1) ~= num_steps, M_m = M_m.'; end
if size(p, 1) ~= num_steps,   p   = p.';   end

% Instantaneous mechanical power per wheel: P_i(t) = |M_m,i(t) * p_i(t)| (Watts)
power_wheels = abs(M_m .* p);
power_t      = sum(power_wheels, 2);

% Trapezoidal integration over time (Joules)
wheel_energy = trapz(t, power_wheels, 1);
total_energy = sum(wheel_energy);

% Optional calculation of paper's mathematical control effort \int ||u||^2 dt
if nargin >= 4 && ~isempty(u)
    if size(u, 1) ~= num_steps, u = u.'; end
    ideal_effort = trapz(t, sum(u.^2, 2));
else
    ideal_effort = NaN;
end
end