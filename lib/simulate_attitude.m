function out = simulate_attitude(c, x0, p, opts)
%SIMULATE_ATTITUDE  Closed loop of a sampled controller c (koopman_mpc) on Eqs. (1)-(2).
%   out = simulate_attitude(c, x0, p, opts), x0 = [w; q]. The torque is held over each
%   sample p.ts and the plant is integrated with ode45.
%   opts fields (all optional):
%     tend      final time, default p.tend
%     jc        true inertia, default p.jc
%     dist      tau_d(t) added to the plant, default none
%     target    add the target torque of Eq. (59), default false
%     noise     add N(0, p.noise_var) to the measured q_v and w, default false
%   out: t, x (7 x nt), tau (3 x nt-1), mu, steptime (controller wall time per sample)

if nargin < 4, opts = struct(); end
d = struct('tend', p.tend, 'jc', p.jc, 'dist', @(t) zeros(3, 1), 'target', false, 'noise', false);
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(opts, f{i}) || isempty(opts.(f{i})), opts.(f{i}) = d.(f{i}); end
end

t = 0:p.ts:opts.tend; nt = numel(t);
x = zeros(7, nt); x(:, 1) = x0;
tau = zeros(3, nt - 1); mu = zeros(1, nt - 1); st = zeros(1, nt - 1);
ode = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
sd = sqrt(p.noise_var);
for k = 1:nt - 1
    xm = x(:, k);
    if opts.noise
        xm(1:3) = xm(1:3) + sd*randn(3, 1);
        xm(5:7) = xm(5:7) + sd*randn(3, 1);
        xm(4:7) = xm(4:7)/norm(xm(4:7));
    end
    t0 = tic;
    [tau(:, k), c, info] = c.step(c, xm);
    st(k) = toc(t0);
    mu(k) = info.mu;
    tk = t(k); uk = tau(:, k);
    [~, xx] = ode45(@(s, y) attitude_ode(y, uk + opts.dist(tk + s) + ...
        opts.target*target_torque(y, p), opts.jc), [0 p.ts/2 p.ts], x(:, k), ode);
    x(:, k + 1) = xx(end, :)';
    x(4:7, k + 1) = x(4:7, k + 1)/norm(x(4:7, k + 1));
end
out = struct('t', t, 'x', x, 'tau', tau, 'mu', mu, 'steptime', st);
end

function ut = target_torque(y, p)
% Eq. (59), u_t = -Kpt q_vt - Kdt w_t, the target rigidly attached so q_t = q and w_t = w [guess]
ut = -p.kpt*y(5:7) - p.kdt*y(1:3);
end
