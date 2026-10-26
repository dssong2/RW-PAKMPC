function ts = settling_time(t, qv, tol)
%SETTLING_TIME  First time after which ||q_v|| stays within tol*||q_v(0)|| (Fig. 4, 5% margin).
%   [guess] The paper says "enter and remain within a 5% error margin" without the reference;
%   the initial attitude error ||q_v(0)|| is used here. Returns NaN if it never settles.
if nargin < 3, tol = 0.05; end
e = vecnorm(qv);
out = find(e > tol*e(1), 1, 'last');
if isempty(out), ts = t(1);
elseif out == numel(t), ts = nan;
else, ts = t(out + 1);
end
end
