function m = rbf_edmd(p, nrbf, seed)
%RBF_EDMD  Data-driven Koopman model by EDMD with radial basis functions, Korda and Mezic [22].
%   Lifted state z = [w; q; psi_1(x); ...; psi_nrbf(x)], size 7 + nrbf (57 for nrbf = 50 as in
%   Sec. V-B), with thin-plate splines psi_i(x) = r^2 log r, r = ||D (x - c_i)||. A and B are fit by
%   least squares to one-step data, z+ = A z + B tau, and C = [I7 0].
%   [guess] The paper does not give the data set. Here: 200 trajectories of 100 samples from
%   random attitudes and |w_i| <= 0.02 rad/s, inputs uniform in [-umax, umax] held per sample,
%   simulated on the true inertia p.jc; centers drawn from the data; w scaled by 1/0.05.
%
%   Same fields as quat_koopman (the input is tau itself, so tau2u and u2tau are identities).

if nargin < 2, nrbf = 50; end
if nargin < 3, seed = 1; end
rs = RandStream('mt19937ar', 'Seed', seed);
ntraj = 200; len = 100;
scale = [ones(3, 1)/0.05; ones(4, 1)];
ode = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

x = zeros(7, ntraj*len); xp = x; u = zeros(3, ntraj*len);
k = 0;
for i = 1:ntraj
    xi = [0.02*(2*rand(rs, 3, 1) - 1); randq(rs)];
    for j = 1:len
        ui = p.umax*(2*rand(rs, 3, 1) - 1);
        [~, xx] = ode45(@(s, y) attitude_ode(y, ui, p.jc), [0 p.ts/2 p.ts], xi, ode);
        k = k + 1;
        x(:, k) = xi; u(:, k) = ui;
        xi = xx(end, :)'; xi(4:7) = xi(4:7)/norm(xi(4:7));
        xp(:, k) = xi;
    end
end
ctr = x(:, randperm(rs, size(x, 2), nrbf)).*scale;

lift = @(xx) [xx; rbf(xx.*scale, ctr)];
z = lift(x); zp = lift(xp);
ab = zp/[z; u];
n = 7 + nrbf;

m.name = 'RBF-based EDMD';
m.n = n; m.nphi = nan; m.ts = p.ts; m.jc0 = p.jc;
m.a = ab(:, 1:n); m.bconst = ab(:, n+1:end);
m.c = [eye(7) zeros(7, nrbf)];
m.t6 = [eye(3) zeros(3, 4); zeros(3, 4) eye(3)];
m.lift = lift;
m.b = @(w, q) m.bconst;
m.tau2u = @(w, tau) tau;
m.u2tau = @(w, u) u;
m.fit = norm(zp - ab*[z; u], 'fro')/norm(zp, 'fro');
end

function psi = rbf(xs, ctr)
% thin-plate spline r^2 log r for each column of xs against each center
r2 = max(sum(xs.^2, 1) + sum(ctr.^2, 1)' - 2*ctr'*xs, 0);
psi = 0.5*r2.*log(r2 + realmin);
end

function q = randq(rs)
q = randn(rs, 4, 1); q = q/norm(q);
if q(1) < 0, q = -q; end
end
