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
%     db        fixed additive Delta B, default zeros
%     adapt     learn Delta B online with Adam, Eqs. (24)-(27) and Algorithm 1, default false
%     accumulate  Algorithm 1 line 11 read as B <- B + Delta B every sample, default false
%     freeze_b  hold the nominal B at B0 = B(q_0, w_0) of the first sample (Algorithm 1 line 2 and
%               Eq. 38) instead of rebuilding B(q_k, w_k) each sample, default false
%
%   With adapt on, each sample the lifted measurement z_k = Phi(x_k) is compared with the nominal
%   one-step prediction z_hat_k = A z_(k-1) + B(x_(k-1)) u_(k-1) (Eq. 23), Dz_k = z_k - z_hat_k,
%   and the last p.nw pairs (u, Dz) form U_k, DZ_k. Once the window is full, one Adam step on the
%   l1-regularized weighted loss (25)-(26) updates Delta B, and the MPC uses B(x_k) + Delta B.
%   [guess] Lambda = diag(lambda^nw, ..., lambda) following the sum in (25); one Adam step per
%   sample; u_(k-1) is the applied (saturated) torque mapped back through Eq. (8).
%
%   [guess] The cost penalizes z - Phi(x_e) rather than z itself, so q -> [1 0 0 0] is the
%   target (Eq. 29a as printed would drive q0 to 0). The input bound |tau_i| <= umax is mapped to
%   u = tau - w^x Jc0 w with w the current measurement, held over the horizon. B is rebuilt each
%   sample from the measured (w, q). The terminal set (28g) is not in the QP form (34), so it is
%   left out, as in the paper's QP.

if nargin < 3, opts = struct(); end
ny = size(m.c, 1);
d = struct('Q', [], 'R', p.R, 'P', [], 'interp', true, 'db', zeros(m.n, 3), 'gamma_e', p.gamma_e, ...
    'adapt', false, 'accumulate', false, 'freeze_b', false);
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

% online Delta B (Algorithm 1)
c.db = opts.db;                             % Adam iterate
c.bacc = zeros(n, 3);                       % accumulated correction when opts.accumulate
c.b0 = [];                                  % B(q_0, w_0) when opts.freeze_b
c.prev = [];                                % lifted state, nominal B and applied u at k - 1
c.uw = zeros(3, 0); c.dzw = zeros(n, 0);    % data window U_k, DZ_k
c.nw = p.nw; c.lam = diag(p.lambda.^(p.nw:-1:1));
c.adam = struct('rho', p.rho, 'gamma', p.gamma, 'b1', p.beta1, 'b2', p.beta2, 'eps', p.eps, ...
    'm', zeros(n, 3), 'v', zeros(n, 3), 't', 0);
c.qpopt = optimoptions('quadprog', 'Display', 'off');
c.step = @step;
end

function [tau, c, info] = step(c, x)
m = c.m; N = c.N; n = m.n;
w = x(1:3); q = x(4:7);
phi = m.lift(x);
if c.opts.freeze_b
    if isempty(c.b0), c.b0 = m.b(w, q); end
    bnom = c.b0;
else
    bnom = m.b(w, q);
end
if c.opts.adapt
    c = update_db(c, phi);
end
if c.opts.accumulate
    c.bacc = c.bacc + c.db;
    b = bnom + c.bacc;
else
    b = bnom + c.db;
end

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
c.prev = struct('z', phi, 'b', bnom + c.opts.accumulate*c.bacc, 'u', m.tau2u(w, tau));
info = struct('mu', mu, 'e', norm(e), 'flag', flag, 'db', norm(c.db));
end

function c = update_db(c, z)
% Eqs. (23)-(27): slide the data window, then one Adam step once it holds nw samples
if isempty(c.prev), return; end
dz = z - (c.m.a*c.prev.z + c.prev.b*c.prev.u);
c.uw = [c.uw(:, max(1, end - c.nw + 2):end) c.prev.u];
c.dzw = [c.dzw(:, max(1, end - c.nw + 2):end) dz];
if size(c.uw, 2) < c.nw, return; end
a = c.adam;
g = -2/c.nw*(c.dzw - c.db*c.uw)*c.lam*c.uw' + a.gamma*sign(c.db);
a.t = a.t + 1;
a.m = a.b1*a.m + (1 - a.b1)*g;
a.v = a.b2*a.v + (1 - a.b2)*g.^2;
mh = a.m/(1 - a.b1^a.t);
vh = a.v/(1 - a.b2^a.t);
c.db = c.db - a.rho*mh./(sqrt(vh) + a.eps);
c.adam = a;
end
