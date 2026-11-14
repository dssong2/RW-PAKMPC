function p = params()
p.Js = diag([405 405 405]); % kg*m2
p.Jt = diag([36.8 37.5 36.8]); % kg*m2
p.ms = 1080; % kg
p.mt = 75; % kg
p.Jc = [986.1 1 -1; 1 573.2 -1; -1 -1 1116]; % kg*m2
p.q0 = [0.5; 0.5; 0.5; 0.5]; % quaternion
p.w0 = [0.02; -0.01; -0.01]; % rad/s
p.q_des = [1;0;0;0]; % quaternion
p.w_des = [0;0;0]; % rad/s
p.umax = 20; % N*m, torque
p.Q = 1e2 * eye(7);
p.R = 1e-2 * eye(3);
p.dt = 0.1; % s
p.sigma = 208.68;
p.gamma_e = 1e6;
p.t_config = 200; % s
end