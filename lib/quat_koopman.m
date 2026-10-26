function m = quat_koopman(nphi, jc0, ts)
%QUAT_KOOPMAN  Physics-based quaternion Koopman model of Sec. III-B, Eqs. (8)-(20).
%   Lifted state z = [w; q; phi_1; ...; phi_nphi], phi_k = Omega(w)^k q, size 4*nphi + 7.
%   Input u = -w^x Jc0 w + tau (Eq. 8). Continuous model z' = Ac z + Bc(w, q) u (Eq. 15),
%   discretized under ZOH as z+ = A z + B u with A = e^(Ac ts), B = G Bc, G = int_0^ts e^(Ac s) ds.
%
%   m.lift(x)      x = [w; q] -> z
%   m.bc(w, q)     Bc(q, w) of Eq. (15)
%   m.b(w, q)      discrete B = G Bc(q, w)
%   m.tau2u(w, tau), m.u2tau(w, u)   Eq. (8) input map

n = 4*nphi + 7;
ac = zeros(n);
for k = 0:nphi - 1
    r = 3 + 4*k;                        % q (k = 0) or phi_k feeds phi_(k+1)
    ac(r + (1:4), r + 4 + (1:4)) = eye(4);
end
e = expm([ac eye(n); zeros(n, 2*n)]*ts);

m.name = 'Quaternion-based';
m.n = n; m.nphi = nphi; m.ts = ts; m.jc0 = jc0; m.jinv = inv(jc0);
m.ac = ac; m.a = e(1:n, 1:n); m.g = e(1:n, n+1:end);
m.c = [eye(7) zeros(7, 4*nphi)];
m.lift = @(x) lift(x, nphi);
m.bc = @(w, q) bc(w, q, m.jinv, nphi, n);
m.b = @(w, q) m.g*bc(w, q, m.jinv, nphi, n);
m.tau2u = @(w, tau) -cross(w, jc0*w) + tau;
m.u2tau = @(w, u) u + cross(w, jc0*w);
end

function z = lift(x, nphi)
w = x(1:3); q = x(4:7);
o = omega_mat(w);
z = [w; q; zeros(4*nphi, 1)];
phi = q;
for k = 1:nphi
    phi = o*phi;
    z(3 + 4*k + (1:4)) = phi;
end
end

function b = bc(w, q, jinv, nphi, n)
% b_k^(c,j) = sum_i Omega^(i-1) Omega(J_j) Omega^(k-i) q, built with
% s_k = Omega s_(k-1) + Omega(J_j) phi_(k-1), s_0 = 0, phi_0 = q
b = zeros(n, 3);
b(1:3, :) = jinv;
o = omega_mat(w);
for j = 1:3
    oj = omega_mat(jinv(:, j));
    s = zeros(4, 1); phi = q;
    for k = 1:nphi
        s = o*s + oj*phi;
        phi = o*phi;
        b(3 + 4*k + (1:4), j) = s;
    end
end
end
