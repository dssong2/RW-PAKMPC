%CHECK_QUAT_KOOPMAN  Step 3 check: dPhi/dt along the true dynamics = Ac Phi + Bc(x) u.
%   The true derivative of the lift is computed with the complex step,
%   dPhi/dt = imag(Phi(x + i h f(x, tau)))/h, which is exact to machine precision because Phi is a
%   polynomial in x. Every block except the last should match the model; the last block should
%   differ by exactly the dropped term phi_(Nphi+1) = Omega(w) phi_Nphi.

p = params(); hp = hyperparams();
nphi = hp.N_phi;
rs = RandStream('mt19937ar', 'Seed', 1);

% test states: the paper's initial state, then random states with faster spins
xs = [p.w0; p.q0];
for i = 1:5
    q = randn(rs, 4, 1);
    xs(:, end + 1) = [0.2*(2*rand(rs, 3, 1) - 1); q/norm(q)];
end
taus = [zeros(3, 1) p.umax*(2*rand(rs, 3, 5) - 1)];

fprintf('--- J_c0 = J_c (model should be exact up to truncation) ---\n');
m = quat_koopman(nphi, p.Jc, p.dt);
names = [{'w', 'q'}, arrayfun(@(k) sprintf('phi_%d', k), 1:nphi, 'UniformOutput', false)];
fprintf('%-6s %-8s', 'state', '|w|');
fprintf(' %9s', names{:}); fprintf(' %12s\n', 'last-phi_N+1');
for i = 1:size(xs, 2)
    x = xs(:, i); w = x(1:3); tau = taus(:, i);
    f = attitude_ode(x, tau, p.Jc);
    dz_true = imag(m.lift(x + 1i*1e-20*f))/1e-20;
    u = m.tau2u(w, tau);
    z = m.lift(x);
    dz_model = m.ac*z + m.bc(x)*u;
    err = dz_true - dz_model;
    fprintf('%-6d %-8.3f', i, norm(w));
    for b = 1:nphi + 2
        r = blk(b);
        fprintf(' %9.1e', norm(err(r))/max(norm(dz_true(r)), eps));
    end
    phinext = omega_mat(w)*z(blk(nphi + 2));
    fprintf(' %12.1e\n', norm(err(blk(nphi + 2)) - phinext));
end
fprintf('(columns are relative errors per block; last column is |residual - phi_(Nphi+1)|)\n');

fprintf('\n--- other checks ---\n');
x = xs(:, 1); w = x(1:3); z = m.lift(x);
fprintf('|phi_k| vs (|w|/2)^k at x0: max rel diff %.1e\n', ...
    max(arrayfun(@(k) abs(norm(z(blk(k + 2))) - (norm(w)/2)^k)/(norm(w)/2)^k, 1:nphi)));
tau = taus(:, 2);
fprintf('u2tau(tau2u(tau)) round trip: %.1e\n', norm(m.u2tau(w, m.tau2u(w, tau)) - tau));
fprintf('size: n = %d, Ac %s, Bc %s, Ac nilpotent: %d\n', m.n, mat2str(size(m.ac)), ...
    mat2str(size(m.bc(x))), nnz(m.ac^(nphi + 1)) == 0);

fprintf('\n--- J_c0 = J_s (Sec. V-C model error) ---\n');
ms = quat_koopman(nphi, p.Js, p.dt);
x = xs(:, 1); w = x(1:3); tau = taus(:, 2);
f = attitude_ode(x, tau, p.Jc);
dz_true = imag(ms.lift(x + 1i*1e-20*f))/1e-20;
dz_model = ms.ac*ms.lift(x) + ms.bc(x)*ms.tau2u(w, tau);
err = dz_true - dz_model;
fprintf('w block: true wdot %s, model wdot %s\n', mat2str(dz_true(1:3)', 3), mat2str(dz_model(1:3)', 3));
fprintf('relative error per block:'); fprintf(' %.1e', arrayfun(@(b) norm(err(blk(b)))/norm(dz_true(blk(b))), 1:nphi + 2)); fprintf('\n');

function r = blk(b)
% rows of block b: 1 = w, 2 = q, 2 + k = phi_k
if b == 1, r = 1:3; else, r = 3 + 4*(b - 2) + (1:4); end
end
