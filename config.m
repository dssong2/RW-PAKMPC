% config.m

% Servicer parameters
ms = 1080;                       % servicer mass (kg)
js = diag([405, 405, 405]);      % servicer inertia tensor (kg*m^2)

% Target parameters (nominal, before parameter sweeps)
mt = 75;                         % target mass (kg)
jt = diag([36.8, 37.5, 36.8]);   %  target inertia tensor (kg*m^2)

% Combined Spacecraft Nominal Post-Capture Inertia (Paper Table I)
jc = [ 986.1,    1.0,   -1.0;
    1.0,  573.2,   -1.0;
    -1.0,   -1.0, 1116.0 ];  % Combined inertia tensor J_c (kg*m^2)

% Initial Tumbling State Post-Capture
q0 = [0.5; 0.5; 0.5; 0.5];       % Initial attitude quaternion [q0; q1; q2; q3] (scalar-first)
w0 = [0.02; -0.01; -0.01];       % Initial body angular rates [wx; wy; wz] (rad/s)

% Desired Stabilized State Reference
q_ref = [1; 0; 0; 0];            % Reference quaternion
w_ref = [0; 0; 0];               % Reference angular rate (rad/s)

% Environmental disturbance & sensor noise toggles
ts        = 0.1;                 % Sampling & discrete solver step Ts (s)
tend      = 200.0;               % Simulation duration (s)
noise_var = 1e-4;                % Sensor noise variance N(0, 1e-4) on w and q_v
noise_on  = 1;                   % 1 = on, 0 = off
dist_on   = 1;                   % 1 = on, 0 = off (unmodeled sinusoidal body disturbance)
target_on = 1;                   % 1 = on, 0 = off (competitive debris reaction torque)

% Target Competitive Controller Gains (Eq. 59)
kpt = 0.02 * jt;                 % Proportional gain matrix
kdt = 0.05 * jt;                 % Derivative gain matrix

% Paper Ideal Actuator Limit
umax = 20.0;                     % Ideal direct 3-axis torque saturation limit per axis (N*m)


%% 2. Tetrahedral Reaction Wheel Actuator Parameters (Sec. IV)
% Tetrahedral Geometry (Pyramid angle beta = 54.74 deg)
beta = acos(1 / sqrt(3));        % beta = 54.7356 deg
wt   = [ sin(beta), -sin(beta),         0,          0;
    cos(beta),  cos(beta),  cos(beta),  cos(beta);
    0,          0,  sin(beta), -sin(beta) ]; % Spin-axis matrix W_T (3x4)
wp   = pinv(wt);                 % Pseudo-inverse steering/allocation matrix (4x3)

% Realistic Actuator Limits & Dynamics (Honeywell HR16-class)
jw       = 0.08;                 % Rotor spin-axis rotational inertia J_w (kg*m^2)
jt_w     = jw / 2;               % Transverse rotor inertia (thin disk approx, kg*m^2)
mw       = 6.0;                  % Mass per reaction wheel rotor (kg)
cth      = 1e-4;                 % Bearing viscous drag c_theta (N*m*s/rad)
taumax   = 0.5;                  % Motor torque saturation limit |tau_w,i| <= tau_max (N*m)
hmax     = 50.0;                 % Rotor momentum saturation limit |h_w,i| <= h_max (N*m*s)

% Actuator Hardware Imperfections (Eq. 5)
sigma_rw = 1e-4 * noise_on;                   % RW torque zero-mean Gaussian noise std (N*m)
b_rw     = 2e-4 * [1, -1, 0.5, -0.5] * noise_on; % RW torque constant bias per wheel (N*m)
sigd_rw  = 1e-5 * noise_on;                   % RW torque random walk drift std (N*m/sqrt(s))


%% 3. Stand-In Controller Gains (Placeholder for PAKMPC)
kp = 40.0;                       % Stand-in attitude proportional gain
kd = 400.0;                      % Stand-in rate damping gain


%% 4. Future PAKMPC Hyperparameters (Table I)
% (Will be used once your team integrates the Koopman MPC block)
N_phi   = 5;                     % Koopman polynomial truncation order
N_w     = 10;                    % Online adaptive sliding window size
N_pred  = 15;                    % MPC prediction horizon
lambda  = 0.965;                 % Exponential forgetting factor
rho     = 1e-7;                  % Adam learning rate
gamma   = 1e-3;                  % Sparsity weight (L1 regularization)
Q_mat   = 1e2 * eye(7);          % State penalty matrix
R_mat   = 1e-2 * eye(3);         % Control effort penalty matrix