function c = incremental_nmpc(p, opts)
%INCREMENTAL_NMPC  INMPC baseline of Sec. V-C, after Chen, Xie and Yu [41]: nonlinear MPC on an
%   incremental model built from measured data.
%   c = incremental_nmpc(p, opts); then each sample [tau, c, info] = c.step(c, x), as koopman_mpc.
%   opts fields (all optional):
%     alpha     low-pass gain on the measured angular acceleration, default 0.5          [guess]
%     maxit     fmincon iterations per sample, default 30                                [guess]
%
%   [paper] "a nonlinear MPC framework that operates directly on an incremental, data-driven model
%   constructed from system data". [guess] Incremental (INDI-type) model: the angular acceleration
%   measured over the last sample, wd_0 = (w_k - w_(k-1))/ts under tau_(k-1), lumps the gyroscopic
%   term, the target torque, the disturbance and the inertia error; over the horizon
%       w_(i+1) = w_i + ts (wd_0 + Jc0^-1 (tau_i - tau_(k-1))),
%       q_(i+1) = exp(Omega(w_i) ts) q_i,
%   with Jc0 = Js the nominal inertia. The cost matches koopman_mpc: p.Q on [w; q] - [0; 1 0 0 0]
%   for i = 1..N-1, p.R on tau, and the LQR terminal weight P6 of terminal_weight on [w; q_v] at
%   i = N. fmincon (SQP) solves for tau_0..tau_(N-1) in [-umax, umax], warm-started by the shifted
%   previous solution.

if nargin < 2, opts = struct(); end
d = struct('alpha', 0.5, 'maxit', 30);
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(opts, f{i}) || isempty(opts.(f{i})), opts.(f{i}) = d.(f{i}); end
end
ms = struct('jc0', p.jc0, 'ts', p.ts, 'c', eye(7), 't6', [eye(3) zeros(3, 4); zeros(3, 4) eye(3)]);
[~, c.p6] = terminal_weight(ms, p);
c.N = p.N; c.ts = p.ts; c.ji = inv(p.jc0); c.umax = p.umax;
c.Q = p.Q; c.R = p.R; c.alpha = opts.alpha;
c.wd = zeros(3, 1); c.wprev = []; c.tauprev = zeros(3, 1);
c.U = zeros(3, p.N);
c.opt = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', ...
    'MaxIterations', opts.maxit, 'SpecifyObjectiveGradient', false);
c.step = @step;
end

function [tau, c, info] = step(c, x)
w = x(1:3); q = x(4:7);
if ~isempty(c.wprev)
    c.wd = (1 - c.alpha)*c.wd + c.alpha*(w - c.wprev)/c.ts;
end
u0 = [c.U(:, 2:end) c.U(:, end)];
lb = -c.umax*ones(3, c.N); ub = -lb;
cost = @(u) rollout(c, w, q, reshape(u, 3, c.N));
[u, ~, flag] = fmincon(cost, u0(:), [], [], [], [], lb(:), ub(:), [], c.opt);
c.U = reshape(u, 3, c.N);
tau = min(max(c.U(:, 1), -c.umax), c.umax);
c.wprev = w; c.tauprev = tau;
info = struct('mu', 0, 'flag', flag, 'db', 0);
end

function j = rollout(c, w, q, u)
qr = [1; 0; 0; 0];
j = 0;
for i = 1:c.N
    du = u(:, i);
    j = j + du'*c.R*du;
    a = norm(w)*c.ts/2;
    if a > 0
        q = cos(a)*q + sin(a)/a*c.ts*omega_mat(w)*q;
    end
    w = w + c.ts*(c.wd + c.ji*(du - c.tauprev));
    if i < c.N
        e = [w; q - qr];
        j = j + e'*c.Q*e;
    else
        e = [w; q(2:4)];
        j = j + e'*c.p6*e;
    end
end
end
