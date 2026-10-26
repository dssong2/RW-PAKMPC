function out = simulate_attitude(c, x0, p, opts)
%SIMULATE_ATTITUDE  Closed loop of a sampled controller c (koopman_mpc) on Eqs. (1)-(2).
%   out = simulate_attitude(c, x0, p, opts), x0 = [w; q]. The torque is held over each
%   sample p.ts and the plant is integrated with ode45.
%   opts fields (all optional):
%     tend      final time, default p.tend
%     jc        true inertia, default p.jc
%     dist      tau_d(t) added to the plant, default none
%     target    add the target torque of Eq. (59), default false
%     target_ref  'initial' (target holds its initial attitude) or 'identity', default 'initial'
%     seed      noise seed, default 1
%     noise     add N(0, p.noise_var) to the measured q_v and w, default false
%   out: t, x (7 x nt), tau (3 x nt-1), mu, db (||Delta B||_2), steptime (controller wall time)

if nargin < 4, opts = struct(); end
d = struct('tend', p.tend, 'jc', p.jc, 'dist', @(t) zeros(3, 1), 'target', false, 'noise', false, ...
    'seed', 1, 'target_ref', 'initial');
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(opts, f{i}) || isempty(opts.(f{i})), opts.(f{i}) = d.(f{i}); end
end

t = 0:p.ts:opts.tend; nt = numel(t);
qt = x0(4:7);
if strcmp(opts.target_ref, 'identity'), qt = [1; 0; 0; 0]; end
x = zeros(7, nt); x(:, 1) = x0;
tau = zeros(3, nt - 1); mu = zeros(1, nt - 1); db = zeros(1, nt - 1); st = zeros(1, nt - 1);
rs = RandStream('mt19937ar', 'Seed', opts.seed);
ode = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
sd = sqrt(p.noise_var);
for k = 1:nt - 1
    xm = x(:, k);
    if opts.noise
        xm(1:3) = xm(1:3) + sd*randn(rs, 3, 1);
        xm(5:7) = xm(5:7) + sd*randn(rs, 3, 1);
        xm(4:7) = xm(4:7)/norm(xm(4:7));
    end
    t0 = tic;
    [tau(:, k), c, info] = c.step(c, xm);
    st(k) = toc(t0);
    mu(k) = info.mu;
    if isfield(info, 'db'), db(k) = info.db; end
    tk = t(k); uk = tau(:, k);
    [~, xx] = ode45(@(s, y) attitude_ode(y, uk + opts.dist(tk + s) + ...
        opts.target*target_torque(y, qt, p), opts.jc), [0 p.ts/2 p.ts], x(:, k), ode);
    x(:, k + 1) = xx(end, :)';
    x(4:7, k + 1) = x(4:7, k + 1)/norm(x(4:7, k + 1));
end
out = struct('t', t, 'x', x, 'tau', tau, 'mu', mu, 'db', db, 'steptime', st);
end

function ut = target_torque(y, q0, p)
% Eq. (59), u_t = -Kpt q_vt - Kdt w_t, q_vt the vector part of q0^* (x) q and w_t = w.
% [guess] As in build_pakmpc_multibody, q0 is the target's initial attitude by default.
q = y(4:7); a = [q0(1); -q0(2:4)];
qv = a(1)*q(2:4) + q(1)*a(2:4) + cross(a(2:4), q(2:4));
ut = -p.kpt*qv - p.kdt*y(1:3);
end
