function m = quat_koopman(nphi, Jc0, ts)
n = 4*nphi + 7;
Ji = inv(Jc0);

% A_c, Eq. (15)
ac = zeros(n);
for k = 0:nphi - 1
    r = 3 + 4*k + (1:4);
    ac(r, r + 4) = eye(4);
end

% ZOH discretization, Eq. (20): expm([Ac I; 0 0] Ts) = [A G; 0 I]
e = expm([ac eye(n); zeros(n, 2*n)]*ts);
g = e(1:n, n+1:2*n);

m.n = n;
m.ts = ts;
m.ac = ac;
m.a = e(1:n, 1:n);                          % A = e^(Ac Ts)
m.g = g;                                    % G = integral of e^(Ac tau), 0 to Ts
m.c = [eye(7) zeros(7, 4*nphi)];            % Cc = [I7 0], Eq. (15); C = Cc, Eq. (20)
m.Ji = Ji;
m.Jc0 = Jc0;
m.lift = @(x) lift(x, nphi);
m.bc   = @(x) bmat(x, nphi, Ji);
m.b    = @(x) g*bmat(x, nphi, Ji);          % B(x) = G Bc(x), Eq. (20)
m.u2tau = @(w, u) u + cross(w, Jc0*w);
m.tau2u = @(w, tau) tau - cross(w, Jc0*w);

% Claude test
A2 = zeros(n); G2 = zeros(n); Ak = eye(n);
for k = 0:nphi + 1
    A2 = A2 + Ak*ts^k/factorial(k);
    G2 = G2 + Ak*ts^(k+1)/factorial(k+1);
    Ak = Ak*m.ac;
end
norm(m.a - A2), norm(m.g - G2)        % both about 1e-16

end

function z = lift(x, nphi)
w = x(1:3); phi = x(4:7);
O = omega_mat(w);
z = zeros(4*nphi + 7, 1); z(1:7) = x;
for k = 1:nphi
    phi = O*phi;
    z(3 + 4*k + (1:4)) = phi;
end
end

function b = bmat(x, nphi, Ji)
w = x(1:3);
q = x(4:7);
O = omega_mat(w);
b = zeros(4*nphi + 7, 3);
b(1:3, :) = Ji;
for j = 1:3
    M = omega_mat(Ji(:, j));             % Omega(J_j), linear in its argument
    phi = q;                             % phi_(k-1), starting at phi_0 = q
    s = M*q;                             % b_1
    b(8:11, j) = s;
    for k = 2:nphi
        phi = O*phi;
        s = O*s + M*phi;
        b(3 + 4*k + (1:4), j) = s;
    end
end
end