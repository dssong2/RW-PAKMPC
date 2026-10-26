function [pz, p6, k6] = terminal_weight(m, p, q6, r)
%TERMINAL_WEIGHT  Terminal weight P on z - z_ref for koopman_mpc.
%   [guess] Replaces the LMI (42)-(43). The lifted pair (A, B0) is not stabilizable (A has every
%   eigenvalue at 1 and rank [A - I, B0] < n), so (40) has no solution with Q > 0 on the lifted
%   state. Instead P6 solves the DARE of the linearization at the target,
%       w+ = w + ts Jc0^-1 u,  q_v+ = q_v + ts/2 w + ts^2/4 Jc0^-1 u,
%   with weights q6 (default p.Q(1) I6) and r (default p.R), and is mapped to the lifted state
%   through the output, P = C' T6' P6 T6 C, so V_F(z_N) is the LQR cost-to-go of the prediction end.

if nargin < 3 || isempty(q6), q6 = p.Q(1)*eye(6); end
if nargin < 4 || isempty(r), r = p.R; end
ji = inv(m.jc0);
ac = [zeros(3, 6); 0.5*eye(3) zeros(3)];
bc = [ji; zeros(3)];
sysd = c2d(ss(ac, bc, eye(6), 0), m.ts);
[p6, k6] = idare(sysd.A, sysd.B, q6, r);
pz = m.c'*m.t6'*p6*m.t6*m.c;
end
