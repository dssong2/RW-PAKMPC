function xd = attitude_ode(x, tau, jc)
%ATTITUDE_ODE  Eqs. (1)-(2) with x = [w; q], tau the total torque (control + disturbance).
w = x(1:3); q = x(4:7);
xd = [jc\(-cross(w, jc*w) + tau); omega_mat(w)*q];
end
