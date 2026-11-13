function [ss_error, t_conv, qn] = calc_convergence(t, q, settle_tol)
% CALC_CONVERGENCE Computes steady-state quaternion vector error and settling time.
%
% Inputs:
%   t          - (N x 1) Simulation time vector (s)
%   q          - (N x 4) Quaternion trajectory [q0, q1, q2, q3]
%   settle_tol - Relative threshold for settling (default: 0.05 for 5%)
%
% Outputs:
%   ss_error   - Final steady-state quaternion vector norm ||q_v(end)||_2
%   t_conv     - Time (s) when ||q_v(t)||_2 permanently enters settle_tol * ||q_v(0)||_2
%   qn         - (N x 1) Time history of ||q_v(t)||_2

if nargin < 3 || isempty(settle_tol)
    settle_tol = 0.05;
end

q_v = q(:, 2:4);
qn  = vecnorm(q_v, 2, 2);

% Steady-state error at the end of the simulation
ss_error = qn(end);

% Last index where attitude error exceeds the settling threshold
threshold = settle_tol * qn(1);
k = find(qn > threshold, 1, 'last');

if isempty(k)
    t_conv = t(1);
elseif k < length(t)
    t_conv = t(k + 1);
else
    t_conv = NaN; % Did not converge within simulation time
end
end