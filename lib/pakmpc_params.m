function p = pakmpc_params()
%PAKMPC_PARAMS  Parameters of Sun et al., Sec. V-A, V-C and Table I.
%   [paper] = quoted from the paper, [guess] = not given in the paper, chosen here

% ---------------- [paper] Sec. V-A ----------------
p.ms = 1080; p.mt = 75;                         % kg
p.js = diag([405 405 405]);                     % servicer, kg m^2
p.jt = diag([36.8 37.5 36.8]);                  % target, kg m^2
p.jc = [986.1 1 -1; 1 573.2 -1; -1 -1 1116];    % combined, kg m^2
p.q0 = [0.5; 0.5; 0.5; 0.5];                    % scalar first
p.w0 = [0.02; -0.01; -0.01];                    % rad/s
p.umax = 20;                                    % |tau_i| <= 20 N m
p.ts = 0.1; p.tend = 200;                       % s
p.Q = 1e2*eye(7); p.R = 1e-2*eye(3);
p.sigma = 208.68;                               % terminal set level from (43)
p.gamma_e = 1e6;                                % (29c)

% ---------------- [paper] Table I ----------------
p.nphi = 5; p.nw = 10; p.N = 15;
p.lambda = 0.965; p.rho = 1e-7; p.gamma = 1e-3;
p.beta1 = 0.9; p.beta2 = 0.999; p.eps = 1e-8;

% ---------------- [paper] Sec. V-C ----------------
p.jc0 = p.js;                                   % nominal inertia for the adaptive runs
p.kpt = 0.02*p.jt; p.kdt = 0.05*p.jt;           % Eq. (59)
p.noise_var = 1e-8;                             % [guess] paper: N(0, 1e-4 I3) on q_v and w. Read as a
                                                %   variance (std 0.01) the torque in Fig. 8 is noise-
                                                %   dominated and ||Delta B|| drifts up instead of
                                                %   decaying as in Fig. 11; std 1e-4 reproduces both
p.dist = @(t) [0.1*sin(0.075*t); -0.2*sin(0.1*t); 0.15*sin(0.1*t)];   % tau_d, N m
end
