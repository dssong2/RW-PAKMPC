function dx = attitude_ode(x, tau, Jc)
w = x(1:3);
q = x(4:7);
dw = Jc \ (-cross(w, Jc*w) + tau);
dq = omega_mat(w) * q;
dx = [dw; dq];
end