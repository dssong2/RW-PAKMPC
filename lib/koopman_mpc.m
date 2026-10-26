function c = koopman_mpc(m, p, opts)
%KOOPMAN_MPC  Koopman MPC with interpolated initial state, Sec. IV-A, Eqs. (28)-(34).
%   c = koopman_mpc(m, p, opts) sets up the controller for model m (quat_koopman, dcm_koopman
%   or rbf_edmd) with parameters p (pakmpc_params). Then, each sample,
%       [tau, c, info] = c.step(c, x)
%   with x = [w; q] the measured state, returns the torque to apply.
%
%   Decision variable v = [U; mu] (Eq. 30), z_0|k = Phi(x_k) + mu e_k (Eq. 31),
%   cost (29) written as the QP (34) with H = 2(Sv' Qb Sv + Su' Rb Su + Hreg), f = 2 Sv' Qb (Sc - Zr).
%   opts fields (all optional):
%     Q, R      output and input weights, default p.Q, p.R (Q is resized for a DCM output)
%     P         terminal weight on z - z_ref, default terminal_weight(m, p, Q6)     [guess]
%     interp    use the interpolated initial state, default true
%     db        additive Delta B, default zeros (non-adaptive)
%
%   [guess] The cost penalizes z - Phi(x_e) rather than z itself, so q -> [1 0 0 0] is the
%   target (Eq. 29a as printed would drive q0 to 0). The input bound |tau_i| <= umax is mapped to
%   u = tau - w^x Jc0 w with w the current measurement, held over the horizon. B is rebuilt each
%   sample from the measured (w, q). The terminal set (28g) is not in the QP form (34), so it is
%   left out, as in the paper's QP.

if nargin < 3, opts = struct(); end
ny = size(m.c, 1);
d = struct('Q', [], 'R', p.R, 'P', [], 'interp', true, 'db', zeros(m.n, 3), 'gamma_e', p.gamma_e);
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(opts, f{i}) || isempty(opts.(f{i})), opts.(f{i}) = d.(f{i}); end
end
if isempty(opts.Q)
    opts.Q = p.Q(1)*eye(ny);                % [guess] same 1e2 weight on every output of a DCM model
end
qz = m.c'*opts.Q*m.c;
if isempty(opts.P), opts.P = terminal_weight(m, p); end

N = p.N; n = m.n;
c.m = m; c.N = N; c.opts = opts; c.umax = p.umax;
c.zref = m.lift([0; 0; 0; 1; 0; 0; 0]);
c.zr = repmat(c.zref, N + 1, 1);

% prediction matrices: Z = Ab z0 + Bb U (Eq. 32)
c.ab = zeros((N + 1)*n, n);
ak = eye(n);
for i = 0:N
    c.ab(i*n + (1:n), :) = ak;
    ak = m.a*ak;
end
c.qb = blkdiag(kron(eye(N), qz), opts.P);
c.rb = kron(eye(N), opts.R);
c.su = [eye(3*N) zeros(3*N, 1)];
c.z1 = [];                                  % z*_1|k-1, for e_k
c.qpopt = optimoptions('quadprog', 'Display', 'off');
c.step = @step;
end

function [tau, c, info] = step(c, x)
m = c.m; N = c.N; n = m.n;
w = x(1:3); q = x(4:7);
phi = m.lift(x);
b = m.b(w, q) + c.opts.db;

% Bb, block lower triangular: block (i, j) = A^(i-j) B for i > j
col = c.ab(1:N*n, :)*b;
bb = zeros((N + 1)*n, 3*N);
for j = 1:N
    bb(j*n + 1:end, 3*(j - 1) + (1:3)) = col(1:(N - j + 1)*n, :);
end

if isempty(c.z1) || ~c.opts.interp
    e = zeros(n, 1);
else
    e = phi - c.z1;
end
sv = [bb c.ab*e];
sc = c.ab*phi;
hreg = zeros(3*N + 1);
hreg(end) = c.opts.gamma_e*(e'*e);
h = 2*(sv'*c.qb*sv + c.su'*c.rb*c.su + hreg);
h = (h + h')/2;
f = 2*sv'*c.qb*(sc - c.zr);

% bounds: |u + w^x Jc0 w| <= umax, 0 <= mu <= 1 (mu fixed at 0 when e = 0)
g = m.u2tau(w, zeros(3, 1));
lb = [repmat(-c.umax - g, N, 1); 0];
ub = [repmat(c.umax - g, N, 1); double(c.opts.interp && any(e))];
[v, ~, flag] = quadprog(h, f, [], [], [], [], lb, ub, [], c.qpopt);
if flag <= 0
    v = [zeros(3*N, 1); 0];
end

u = v(1:3);
mu = v(end);
z0 = phi + mu*e;
c.z1 = m.a*z0 + b*u;
tau = min(max(m.u2tau(w, u), -c.umax), c.umax);
info = struct('mu', mu, 'e', norm(e), 'flag', flag);
end
