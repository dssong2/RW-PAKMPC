function m = dcm_koopman(nphi, jc0, ts)
%DCM_KOOPMAN  DCM-based Koopman model of Chen, Shan and Wen [18], built the same way as
%   quat_koopman from R' = R w^x: phi_k = vec(R (w^x)^k), with
%   phi_k' = phi_(k+1) + sum_i vec(R (w^x)^(i-1) (J u)^x (w^x)^(k-i)).
%   Lifted state z = [w; vec(R); phi_1; ...; phi_nphi], size 9*nphi + 12.
%   [guess] The paper gives 9(nphi + 1) in the text and 9*nphi in Table II; [18] is not
%   reproduced in the paper, so w is kept in the state here as in the quaternion model.
%
%   Same fields as quat_koopman; m.lift takes x = [w; q] and m.b/m.bc take (w, q).

n = 9*nphi + 12;
ac = zeros(n);
for k = 0:nphi - 1
    r = 3 + 9*k;                        % vec(R) (k = 0) or phi_k feeds phi_(k+1)
    ac(r + (1:9), r + 9 + (1:9)) = eye(9);
end
e = expm([ac eye(n); zeros(n, 2*n)]*ts);

m.name = 'DCM-based';
m.n = n; m.nphi = nphi; m.ts = ts; m.jc0 = jc0; m.jinv = inv(jc0);
m.ac = ac; m.a = e(1:n, 1:n); m.g = e(1:n, n+1:end);
m.c = [eye(12) zeros(12, 9*nphi)];
m.lift = @(x) lift(x, nphi);
m.bc = @(w, q) bc(w, q2r(q), m.jinv, nphi, n);
m.b = @(w, q) m.g*bc(w, q2r(q), m.jinv, nphi, n);
m.tau2u = @(w, tau) -cross(w, jc0*w) + tau;
m.u2tau = @(w, u) u + cross(w, jc0*w);
end

function z = lift(x, nphi)
w = x(1:3); r = q2r(x(4:7));
wx = skew(w);
z = [w; r(:); zeros(9*nphi, 1)];
phi = r;
for k = 1:nphi
    phi = phi*wx;
    z(3 + 9*k + (1:9)) = phi(:);
end
end

function b = bc(w, r, jinv, nphi, n)
% T_k = sum_i W^(i-1) S W^(k-i), built with T_k = T_(k-1) W + W^(k-1) S, T_0 = 0
b = zeros(n, 3);
b(1:3, :) = jinv;
wx = skew(w);
for j = 1:3
    sj = skew(jinv(:, j));
    t = zeros(3); wp = eye(3);
    for k = 1:nphi
        t = t*wx + wp*sj;
        wp = wp*wx;
        s = r*t;
        b(3 + 9*k + (1:9), j) = s(:);
    end
end
end
