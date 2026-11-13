function [peak_ratio, max_hw_per_wheel, is_saturated, sat_time] = check_momentum_limits(t, hw, hmax)
% CHECK_MOMENTUM_LIMITS Evaluates reaction wheel momentum management and saturation.
%
% Inputs:
%   t    - (N x 1) Simulation time vector (s)
%   hw   - (N x 4) Reaction wheel angular momenta h_{w,i} = J_w * p_i (N*m*s)
%   hmax - Scalar maximum momentum capacity per wheel (N*m*s)
%
% Outputs:
%   peak_ratio       - Max |h_w,i| / h_max across all wheels (dimensionless ratio)
%   max_hw_per_wheel - (1 x 4) Peak absolute momentum per wheel (N*m*s)
%   is_saturated     - Boolean true if any wheel reaches >= 99.9% of h_max
%   sat_time         - Total duration (s) spent in momentum saturation

t = t(:);
if size(hw, 1) ~= length(t), hw = hw.'; end

abs_hw           = abs(hw);
max_hw_per_wheel = max(abs_hw, [], 1);
peak_ratio       = max(max_hw_per_wheel) / hmax;

% Check saturation condition (within 0.1% of hmax)
sat_mask     = any(abs_hw >= 0.999 * hmax, 2);
is_saturated = any(sat_mask);

if is_saturated && length(t) > 1
    dt       = [diff(t); 0];
    sat_time = sum(dt(sat_mask));
else
    sat_time = 0.0;
end
end