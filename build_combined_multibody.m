%BUILD_COMBINED_MULTIBODY  Unified Simscape Multibody model of the combined spacecraft
%   from Sun et al. (PAKMPC paper, Sec. V-A & Fig. 5) integrated with a tetrahedral
%   4-wheel reaction wheel array inside a translucent servicer body.
%
%   Features:
%   1) Exact post-capture inertia tensor Jc, servicer/target bodies, arms, panels, ring.
%   2) Translucent servicer body (Opacity = 0.35) so internal reaction wheels are visible.
%   3) Tetrahedral 4-wheel reaction wheel assembly inside the servicer with red off-axis
%      spoke markers on each rotor so wheel spinning is clearly visible in 3D.
%   4) Runs three automatic simulations:
%      - Run 0: Free-drift total angular momentum conservation check (H_body + H_wheels).
%      - Run 1: Idealized Direct-Torque Actuation (paper's black-box baseline).
%      - Run 2: Realistic Reaction Wheel Actuation (with rotor inertia, friction, saturation,
%               noise, bias, and drift).
%   5) Generates side-by-side comparison plots and RW telemetry plots.

clear; clc;

% ---------------- Toggles & Configuration ----------------
noise_on    = 1;        % Sec. V-C sensor noise + wheel torque noise
dist_on     = 1;        % Sec. V-C external disturbance torque
target_on   = 1;        % Eq. (59) competitive target PD torque
drift_check = true;     % Run free-drift momentum check first

% Reaction wheel sizing preset:
%   'realistic_hr16' -> Realistic Honeywell HR16-class wheels (tau_max = 0.5 N*m, h_max = 50 N*m*s)
%                       Shows clear actuator saturation & slower settling vs paper's 20 N*m ideal!
%   'high_torque'    -> Over-sized wheels matched to paper's 20 N*m authority (tau_max = 15 N*m)
rw_preset = 'realistic_hr16';

v = ver;
if ~any(strcmp({v.Name}, 'Simscape')) || ~any(strcmp({v.Name}, 'Simscape Multibody'))
    error('Need Simscape and Simscape Multibody: Home > Add-Ons > Get Add-Ons.');
end

here = fileparts(mfilename('fullpath'));
if isempty(here) || contains(here, 'Editor_'), here = pwd; end

% ---------------- [paper] Sec. V-A and V-C Spacecraft Parameters ----------------
ms = 1080; mt = 75;                             % Masses (kg)
js = diag([405 405 405]);                       % Servicer inertia (kg m^2)
jt = diag([36.8 37.5 36.8]);                    % Target inertia (kg m^2)
jc = [986.1 1 -1; 1 573.2 -1; -1 -1 1116];      % Combined inertia Jc (kg m^2)
q0 = [0.5; 0.5; 0.5; 0.5];                      % Initial quaternion (scalar first)
w0 = [0.02; -0.01; -0.01];                      % Initial angular velocity (rad/s)
umax = 20;                                      % Paper's idealized torque limit |u_i| <= 20 N m
ts = 0.1; tend = 200;                           % Sampling time and duration (s)
kpt = 0.02*jt; kdt = 0.05*jt;                   % Eq. (59) target competitive PD gains
noise_var = 1e-4;                               % N(0, 1e-4) sensor noise variance
settle_tol = 0.05;                              % 5% settling error margin

% ---------------- Reaction Wheel Parameters (Tetrahedral 4-Wheel Array) ----------------
beta = acos(1/sqrt(3));                         % Tetrahedral angle, 54.74 deg
wt = [sin(beta) -sin(beta)  0          0;
      cos(beta)  cos(beta)  cos(beta)  cos(beta);
      0          0          sin(beta) -sin(beta)]; % Spin axes matrix W_T (3x4)
wp = pinv(wt);                                  % Pseudo-inverse allocation matrix (4x3)

switch rw_preset
    case 'realistic_hr16'
        jw     = 0.08;      % Rotor axial inertia J_w (kg m^2) [HR16-class for ~1155 kg sat]
        mw     = 6.0;       % Rotor mass (kg)
        taumax = 0.5;       % Max motor torque per wheel |M_m,i| <= 0.5 N m (vs 20 N m ideal)
        hmax   = 50;        % Max momentum capacity per wheel |h_w,i| <= 50 N m s
    case 'high_torque'
        jw     = 0.15;      % Larger rotor inertia (kg m^2)
        mw     = 8.0;       % Rotor mass (kg)
        taumax = 15.0;      % High motor torque per wheel (N m)
        hmax   = 100;       % Momentum capacity (N m s)
    otherwise
        error('Unknown rw_preset.');
end
jt_w  = jw / 2;                                 % Rotor transverse inertia (thin disk approx, kg m^2)
cth   = 1e-4;                                   % Bearing viscous drag c_theta (N m s/rad)
sigma_rw = 1e-4 * noise_on;                     % Wheel torque noise std (N m)
b_rw     = 2e-4 * [1 -1 0.5 -0.5] * noise_on;   % Constant wheel torque bias (N m)
sigd_rw  = 1e-5 * noise_on;                     % Wheel torque random-walk drift (N m / sqrt(s))

% ---------------- Spacecraft Geometry & Inertia Accounting ----------------
ls = sqrt(6*js(1,1)/ms);                        % Servicer cube side -> 1.5 m
lt = sqrt(6/mt*(trace(jt) - 2*diag(jt)'));      % Target box -> 1.73 x 1.70 x 1.73 m
gap = 0.2;                                      % Docking ring gap between faces (m)
r = [0; ls/2 + gap + lt(2)/2; 0];               % Target CoM relative to Servicer CoM (m)
c = mt*r/(ms + mt);                             % Combined CoM in Servicer frame (m)
kp = 40; kd = 400;                              % Stand-in PD gains

% Wheel placement inside the 1.5 m servicer cube (tetrahedral cluster around Servicer center)
rw_cluster_rad = 0.32;                          % Radial offset from servicer center (m)
rw_cyl_rad     = 0.16;                          % Visual wheel cylinder radius (m)
rw_cyl_len     = 0.06;                          % Visual wheel cylinder thickness (m)

% Calculate inertia contributed by the 4 reaction wheels about the combined CoM
% so we can subtract the non-axial wheel inertia from the residual appendage inertia
% and keep total combined inertia equal to the paper's Jc!
pa = @(m, d) m*((d'*d)*eye(3) - d*(d'));
j_wheels_locked = zeros(3,3);
j_wheels_trans  = zeros(3,3);
m_wheels_total  = 4 * mw;
for i = 1:4
    pos_w_sv = rw_cluster_rad * wt(:,i);        % Wheel i position in servicer frame
    pos_w_cm = -c + pos_w_sv;                   % Wheel i position in combined CoM frame
    iw_rot_locked = jt_w*(eye(3) - wt(:,i)*wt(:,i)') + jw*(wt(:,i)*wt(:,i)');
    iw_rot_trans  = jt_w*(eye(3) - wt(:,i)*wt(:,i)');
    j_wheels_locked = j_wheels_locked + iw_rot_locked + pa(mw, pos_w_cm);
    j_wheels_trans  = j_wheels_trans  + iw_rot_trans  + pa(mw, pos_w_cm);
end

% Adjust servicer mass so total mass (servicer + wheels + target) remains ms + mt
ms_hull = ms - m_wheels_total;
jst = js + pa(ms_hull, -c) + jt + pa(mt, r - c) + j_wheels_locked;
jres = jc - jst;
% Effective body inertia when wheels are free to spin around their axes:
j_body_eff = jc - (j_wheels_locked - j_wheels_trans);

e = sort(eig((jres + jres')/2));
fprintf('Appendage residual inertia principal moments: [%.2f %.2f %.2f] kg m^2\n', e);
if e(1) <= 0
    error('Residual inertia is not positive definite.');
end

% ---------------- Build Simscape Multibody Model ----------------
mdl = 'pakmpc_rw_combined';
if bdIsLoaded(mdl), close_system(mdl, 0); end
load_system('sm_lib'); load_system('nesl_utility'); load_system('simulink');
new_system(mdl); open_system(mdl);
set_param(mdl, 'Solver', 'ode45', 'RelTol', '1e-8', 'AbsTol', '1e-10', ...
    'MaxStep', '0.01', 'StopTime', num2str(tend));
try, set_param(mdl, 'SimMechanicsOpenEditorOnUpdate', 'on'); catch, end
m = @(x) mat2str(x, 17);
pos = @(x, y) [x y x+50 y+40];

world = add_block('sm_lib/Frames and Transforms/World Frame', [mdl '/World'], 'Position', pos(30, 300));
mech = add_block('sm_lib/Utilities/Mechanism Configuration', [mdl '/Mech'], 'Position', pos(30, 380));
solver = add_block('nesl_utility/Solver Configuration', [mdl '/Solver'], 'Position', pos(30, 460));
if ~sp(mech, {'UniformGravity'}, {'None'}, true)
    sp(mech, {'GravityVector'}, {'[0 0 0]'});
end

six = add_block('sm_lib/Joints/6-DOF Joint', [mdl '/Combined Body Joint'], 'Position', pos(150, 300));
f = fieldnames(get_param(six, 'DialogParameters'));
sp(six, {'SphPositionTargetSpecify'}, {'on'});
sp(six, {'SphPositionTargetRotationMethod'}, {'Quaternion'});
sp(six, {'SphPositionTargetQuaternion'}, {m(q0')});
sp(six, {'SphVelocityTargetSpecify'}, {'on'});
sp(six, {'SphVelocityTargetValueUnits'}, {'rad/s'});
sp(six, {'SphVelocityTargetValue'}, {m(w0')});
sp(six, {'SphVelocityTargetInFollowerFrame'}, {'on'});
sp(six, {'SphTorqueActuationMode'}, {'InputTorque', 'Provided by Input', 'ProvidedByInput'});
sp(six, {'SphActuateTorque'}, {'on'});
sp(six, {'SphActuationFrame'}, {'Follower', 'FollowerFrame'});
sp(six, {'SphSensePosition'}, {'on'});
sp(six, {'SphSenseVelocity'}, {'on'});
sp(six, {'SphSensingFrame'}, {'Follower', 'FollowerFrame', 'SensingFrameFollower'});

ph = get_param(six, 'PortHandles');
if numel(ph.LConn) ~= 2 || numel(ph.RConn) ~= 3
    error('Combined Body Joint ports mismatch.');
end
add_line(mdl, one(world), ph.LConn(1));
add_line(mdl, one(world), one(mech));
add_line(mdl, one(world), one(solver));
cm = ph.RConn(1);           % Combined CoM frame

% ---------------- Spacecraft Bodies & Translucent Servicer ----------------
addinertia(mdl, cm, 'Appendages (arms, panels, ring)', 1e-6, jres, [280 40]);
sv = addrt(mdl, cm, 'Servicer Frame', -c, [], [280 120]);
addinertia(mdl, sv, 'Servicer Hull Inertia', ms_hull, js, [380 120]);

% Translucent Servicer Box (Opacity = 0.35 so internal wheels are clearly visible!)
addbox_translucent(mdl, sv, 'Servicer Body (Translucent)', [0 0 0], [ls ls ls], [0.25 0.55 0.9], 0.35, [460 120]);

tg = addrt(mdl, sv, 'Target Frame', r, [], [380 200]);
addinertia(mdl, tg, 'Target', mt, jt, [480 200]);
addbox_translucent(mdl, tg, 'Target Body', [0 0 0], lt, [0.85 0.7 0.25], 1.0, [560 200]);

% Docking ring between faces
addcyl(mdl, sv, 'Docking Ring', [0 ls/2 0], [0 ls/2 + gap 0], 0.35, [0.35 0.35 0.35], [460 280]);

% Two three-link robotic arms
k = 0;
for s = [-1 1]
    pts = [s*0.3*ls,         ls/2,        0;
           s*(lt(1)/2+0.35), ls/2 + 0.25, 0;
           s*(lt(1)/2+0.45), r(2) - 0.2,  0;
           s*lt(1)/2,        r(2),        0];
    for j = 1:3
        k = k + 1;
        addcyl(mdl, sv, sprintf('Arm %d Link %d', (s+3)/2, j), pts(j,:), pts(j+1,:), 0.05, ...
            [0.3 0.3 0.35], [460 + 110*j, 280 + 55*k]);
    end
end

% Solar panels on Target's +-x sides
for s = [-1 1]
    i = (s+3)/2;
    addcyl(mdl, tg, sprintf('Panel Boom %d', i), [s*lt(1)/2 0 0], [s*(lt(1)/2 + 0.25) 0 0], 0.03, ...
        [0.5 0.5 0.5], [700, 120 + 60*i]);
    addbox_translucent(mdl, tg, sprintf('Solar Panel %d', i), [s*(lt(1)/2 + 0.25 + 0.9) 0 0], [1.8 0.9 0.03], ...
        [0.85 0.65 0.15], 1.0, [820, 120 + 60*i]);
end

% ---------------- Tetrahedral 4-Wheel Assembly Inside Servicer ----------------
mux_p = add_block('simulink/Signal Routing/Mux', [mdl '/pMux'], 'Inputs', '4', 'Position', [820 920 830 1120]);
dmx_m = add_block('simulink/Signal Routing/Demux', [mdl '/MmDemux'], 'Outputs', '4', 'Position', [1320 920 1330 1120]);

for i = 1:4
    y = 900 + 110*(i-1);
    rot_w = [null(wt(:,i)') wt(:,i)];
    if det(rot_w) < 0, rot_w(:,1) = -rot_w(:,1); end
    pos_w = rw_cluster_rad * wt(:,i);

    % Rigid Transform from Servicer Center to Wheel i location & spin axis (+Z)
    rt = add_block('sm_lib/Frames and Transforms/Rigid Transform', sprintf('%s/RW_RT%d', mdl, i), 'Position', pos(280, y));
    sp(rt, {'RotationMethod'}, {'RotationMatrix', 'Rotation Matrix'});
    sp(rt, {'RotationMatrix'}, {m(rot_w)});
    sp(rt, {'TranslationMethod'}, {'Cartesian'});
    sp(rt, {'TranslationCartesianOffset'}, {m(pos_w')});

    % Revolute Joint for Wheel i
    rev = add_block('sm_lib/Joints/Revolute Joint', sprintf('%s/Wheel Joint %d', mdl, i), 'Position', pos(400, y));
    sp(rev, {'TorqueActuationMode'}, {'InputTorque', 'Provided by Input'});
    sp(rev, {'MotionActuationMode'}, {'ComputedMotion', 'Automatically Computed'});
    sp(rev, {'SenseVelocity'}, {'on'});
    if sp(rev, {'DampingCoefficientUnits'}, {'N*m/(rad/s)'}, true)
        sp(rev, {'DampingCoefficient'}, {num2str(cth, 17)});
    else
        sp(rev, {'DampingCoefficient'}, {num2str(cth*pi/180, 17)});
    end

    % Physical Wheel Inertia
    whl = add_block('sm_lib/Body Elements/Inertia', sprintf('%s/Wheel %d', mdl, i), 'Position', pos(540, y));
    sp(whl, {'InertiaType'}, {'Custom', 'CustomInertia'});
    sp(whl, {'Mass'}, {num2str(mw, 17)});
    sp(whl, {'MomentsOfInertia'}, {m([jt_w jt_w jw])});

    % Visual Wheel Cylinder (Dark Gold/Orange Rotor)
    cyl = add_block('sm_lib/Body Elements/Cylindrical Solid', sprintf('%s/Wheel Graphic %d', mdl, i), 'Position', pos(540, y-50));
    sp(cyl, {'CylinderRadius'}, {num2str(rw_cyl_rad)});
    sp(cyl, {'CylinderLength'}, {num2str(rw_cyl_len)});
    sp(cyl, {'Density'}, {'1e-9'});
    sp(cyl, {'GraphicDiffuseColor'}, {'[0.9 0.45 0.1]'}, true);

    % Visual Spin Marker (Bright Red Spoke/Dot offset along X on the rotor face so spin is visible!)
    mkt = add_block('sm_lib/Frames and Transforms/Rigid Transform', sprintf('%s/Marker Offset %d', mdl, i), 'Position', pos(630, y-50));
    sp(mkt, {'TranslationMethod'}, {'Cartesian'});
    sp(mkt, {'TranslationCartesianOffset'}, {m([0.55*rw_cyl_rad, 0, 0])});
    mkb = add_block('sm_lib/Body Elements/Brick Solid', sprintf('%s/Spin Marker %d', mdl, i), 'Position', pos(710, y-50));
    sp(mkb, {'BrickDimensions'}, {m([0.7*rw_cyl_rad, 0.22*rw_cyl_rad, 1.35*rw_cyl_len])});
    sp(mkb, {'Density'}, {'1e-9'});
    sp(mkb, {'GraphicDiffuseColor'}, {'[0.95 0.1 0.1]'}, true);

    pso = add_block('nesl_utility/PS-Simulink Converter', sprintf('%s/p%d', mdl, i), 'Position', pos(680, y));
    sp(pso, {'Unit'}, {'rad/s'});
    spt_w = add_block('nesl_utility/Simulink-PS Converter', sprintf('%s/Mm%d', mdl, i), 'Position', pos(1380, y));
    sp(spt_w, {'Unit'}, {'N*m'});

    pr = get_param(rev, 'PortHandles');
    add_line(mdl, sv, lc(rt, 1));
    add_line(mdl, rc(rt, 1), pr.LConn(1));
    add_line(mdl, pr.RConn(1), one(whl));
    add_line(mdl, pr.RConn(1), one(cyl));
    add_line(mdl, pr.RConn(1), lc(mkt, 1));
    add_line(mdl, rc(mkt, 1), one(mkb));

    add_line(mdl, pr.RConn(2), lc(pso, 1));
    add_line(mdl, op(pso, 1), ip(mux_p, i));
    add_line(mdl, op(dmx_m, i), ip(spt_w, 1));
    add_line(mdl, rc(spt_w, 1), pr.LConn(2));
end

% ---------------- Sensing & Controller ----------------
psq = add_block('nesl_utility/PS-Simulink Converter', [mdl '/Q'], 'Position', pos(280, 520));
psw = add_block('nesl_utility/PS-Simulink Converter', [mdl '/w'], 'Position', pos(280, 590));
sp(psq, {'Unit'}, {'1'});
sp(psw, {'Unit'}, {'rad/s'});
add_line(mdl, ph.RConn(2), lc(psq, 1));
add_line(mdl, ph.RConn(3), lc(psw, 1));

zq = add_block('simulink/Discrete/Zero-Order Hold', [mdl '/ZOHq'], 'SampleTime', num2str(ts), 'Position', pos(400, 520));
zw = add_block('simulink/Discrete/Zero-Order Hold', [mdl '/ZOHw'], 'SampleTime', num2str(ts), 'Position', pos(400, 590));
zp = add_block('simulink/Discrete/Zero-Order Hold', [mdl '/ZOHp'], 'SampleTime', num2str(ts), 'Position', pos(870, 1000));

nz = add_block('simulink/Sources/Random Number', [mdl '/Noise'], 'Mean', m(zeros(1,6)), ...
    'Variance', m(noise_var*ones(1,6)), 'Seed', m(1:6), 'SampleTime', num2str(ts), 'Position', pos(400, 660));
gn = add_block('simulink/Math Operations/Gain', [mdl '/Noise On'], 'Gain', num2str(noise_on), 'Position', pos(480, 660));
sen = add_block('simulink/User-Defined Functions/MATLAB Function', [mdl '/Sensor'], 'Position', [560 520 680 680]);
setscript([mdl '/Sensor'], {
    'function [wm, qm] = sensor(q, w, n)'
    'q = q(:); w = w(:); n = n(:);'
    'wm = w + n(1:3);'
    'qm = [q(1); q(2:4) + n(4:6)];'
    'qm = qm/norm(qm);'});

ctl = add_block('simulink/User-Defined Functions/MATLAB Function', [mdl '/Controller'], 'Position', [740 520 880 640]);
setscript([mdl '/Controller'], {
    'function tau = pd(wm, qm)'
    '% Stand-in controller for PAKMPC u*_0|k'
    sprintf('kp = %.17g; kd = %.17g;', kp, kd)
    's = 1; if qm(1) < 0, s = -1; end'
    'tau = -kp*s*qm(2:4) - kd*wm;'});
sat = add_block('simulink/Discontinuities/Saturation', [mdl '/Torque Limit'], 'UpperLimit', num2str(umax), ...
    'LowerLimit', num2str(-umax), 'Position', pos(920, 540));

% Mode Gains to switch cleanly between Ideal Direct Torque (Paper) and Reaction Wheels
g_ideal = add_block('simulink/Math Operations/Gain', [mdl '/Ideal Direct Mode'], 'Gain', '0', 'Position', pos(1010, 540));
g_rw    = add_block('simulink/Math Operations/Gain', [mdl '/RW Mode'], 'Gain', '1', 'Position', pos(1010, 940));

% ---------------- Reaction Wheel Disturbances & Allocation ----------------
nz_w = add_block('simulink/Sources/Random Number', [mdl '/RW Noise'], 'Mean', m(zeros(1,4)), ...
    'Variance', m((sigma_rw^2 + 1e-30)*ones(1,4)), 'Seed', m(21:24), 'SampleTime', num2str(ts), 'Position', pos(800, 780));
bs_w = add_block('simulink/Sources/Constant', [mdl '/RW Bias'], 'Value', m(b_rw), 'Position', pos(800, 840));
dn_w = add_block('simulink/Sources/Random Number', [mdl '/RW Drift Inc'], 'Mean', m(zeros(1,4)), ...
    'Variance', m((sigd_rw^2*ts + 1e-30)*ones(1,4)), 'Seed', m(31:34), 'SampleTime', num2str(ts), 'Position', pos(720, 900));
di_w = add_block('simulink/Discrete/Discrete-Time Integrator', [mdl '/RW Drift'], 'SampleTime', num2str(ts), ...
    'Position', pos(800, 900));
sp(di_w, {'gainval'}, {num2str(1/ts, 17)});
ds_w = add_block('simulink/Math Operations/Sum', [mdl '/Wheel Disturbance'], 'Inputs', '+++', 'Position', pos(910, 840));
g_wdist = add_block('simulink/Math Operations/Gain', [mdl '/RW Dist On'], 'Gain', num2str(noise_on), 'Position', pos(970, 840));

alloc = add_block('simulink/User-Defined Functions/MATLAB Function', [mdl '/RW Allocation'], 'Position', [1120 920 1270 1040]);
setscript([mdl '/RW Allocation'], {
    'function mm = alloc(u, dw, p)'
    '% Maps desired 3-axis torque u into 4 tetrahedral wheels: M_m = -W_T^+ u - dw'
    'u = u(:); dw = dw(:); p = p(:);'
    ['wp = ' m(wp) ';']
    sprintf('jw = %.17g; taumax = %.17g; hmax = %.17g;', jw, taumax, hmax)
    'mm = -wp*u - dw;'
    'mm = min(max(mm, -taumax), taumax);'
    'mm(abs(jw*p) >= hmax & sign(mm) == sign(p)) = 0;'});

% ---------------- Target Torque, Eq. (59), and External Disturbance ----------------
tgt = add_block('simulink/User-Defined Functions/MATLAB Function', [mdl '/Target Torque'], 'Position', [560 720 680 780]);
setscript([mdl '/Target Torque'], {
    'function ut = target(q, w)'
    'q = q(:); w = w(:);'
    ['a = ' m([q0(1); -q0(2:4)]) ';']
    ['kpt = ' m(kpt) '; kdt = ' m(kdt) ';']
    'qv = a(1)*q(2:4) + q(1)*a(2:4) + cross(a(2:4), q(2:4));'
    'ut = -kpt*qv - kdt*w;'});
gt = add_block('simulink/Math Operations/Gain', [mdl '/Target On'], 'Gain', num2str(target_on), 'Position', pos(720, 730));

clk = add_block('simulink/Sources/Clock', [mdl '/Clock'], 'Position', pos(480, 820));
dst = add_block('simulink/User-Defined Functions/MATLAB Function', [mdl '/Disturbance'], 'Position', [560 810 680 860]);
setscript([mdl '/Disturbance'], {
    'function td = dist(t)'
    'td = [0.1*sin(0.075*t); -0.2*sin(0.1*t); 0.15*sin(0.1*t)];'});
gd = add_block('simulink/Math Operations/Gain', [mdl '/Disturbance On'], 'Gain', num2str(dist_on), 'Position', pos(720, 820));

tsum = add_block('simulink/Math Operations/Sum', [mdl '/External Torque Sum'], 'Inputs', '+++', 'Position', pos(1100, 640));
spt_body = add_block('nesl_utility/Simulink-PS Converter', [mdl '/Body Torque'], 'Position', pos(1180, 640));
sp(spt_body, {'Unit'}, {'N*m'});

% ---------------- Wire Signal Connections ----------------
add_line(mdl, op(psq, 1), ip(zq, 1));
add_line(mdl, op(psw, 1), ip(zw, 1));
add_line(mdl, op(mux_p, 1), ip(zp, 1));
add_line(mdl, op(nz, 1), ip(gn, 1));
add_line(mdl, op(zq, 1), ip(sen, 1));
add_line(mdl, op(zw, 1), ip(sen, 2));
add_line(mdl, op(gn, 1), ip(sen, 3));
add_line(mdl, op(sen, 1), ip(ctl, 1));
add_line(mdl, op(sen, 2), ip(ctl, 2));
add_line(mdl, op(ctl, 1), ip(sat, 1));

% Branch saturated controller output to Ideal Direct path and RW Allocation path
add_line(mdl, op(sat, 1), ip(g_ideal, 1));
add_line(mdl, op(sat, 1), ip(g_rw, 1));

% Wire RW Disturbance & Allocation
add_line(mdl, op(nz_w, 1), ip(ds_w, 1));
add_line(mdl, op(bs_w, 1), ip(ds_w, 2));
add_line(mdl, op(dn_w, 1), ip(di_w, 1));
add_line(mdl, op(di_w, 1), ip(ds_w, 3));
add_line(mdl, op(ds_w, 1), ip(g_wdist, 1));
add_line(mdl, op(g_rw, 1), ip(alloc, 1));
add_line(mdl, op(g_wdist, 1), ip(alloc, 2));
add_line(mdl, op(zp, 1), ip(alloc, 3));
add_line(mdl, op(alloc, 1), ip(dmx_m, 1));

% Wire External Body Torques (Ideal Direct + Target + Disturbance)
add_line(mdl, op(psq, 1), ip(tgt, 1));
add_line(mdl, op(psw, 1), ip(tgt, 2));
add_line(mdl, op(tgt, 1), ip(gt, 1));
add_line(mdl, op(clk, 1), ip(dst, 1));
add_line(mdl, op(dst, 1), ip(gd, 1));
add_line(mdl, op(g_ideal, 1), ip(tsum, 1));
add_line(mdl, op(gt, 1), ip(tsum, 2));
add_line(mdl, op(gd, 1), ip(tsum, 3));
add_line(mdl, op(tsum, 1), ip(spt_body, 1));
add_line(mdl, rc(spt_body, 1), ph.LConn(2));

% ---------------- Workspace Loggers ----------------
logs = {'q_log', psq, 1; 'w_log', psw, 1; 'p_log', mux_p, 1; 'mm_log', alloc, 1; 'u_log', sat, 1};
for i = 1:size(logs, 1)
    lg = add_block('simulink/Sinks/To Workspace', [mdl '/' logs{i,1}], 'VariableName', logs{i,1}, ...
        'SaveFormat', 'Structure With Time', 'Position', pos(1460, 40 + 70*i));
    add_line(mdl, op(logs{i,2}, logs{i,3}), ip(lg, 1));
end

save_system(mdl, fullfile(here, [mdl '.slx']));

% =========================================================================
% RUN 0: FREE-DRIFT ANGULAR MOMENTUM CONSERVATION SANITY CHECK
% =========================================================================
if drift_check
    fprintf('\n--- Run 0: Free-Drift Angular Momentum Conservation Check ---\n');
    set_param(g_ideal, 'Gain', '0'); set_param(g_rw, 'Gain', '0'); set_param(g_wdist, 'Gain', '0');
    set_param(gt, 'Gain', '0'); set_param(gd, 'Gain', '0');
    out0 = sim(mdl, 'StopTime', '60');
    t0 = out0.q_log.time; q_d = sig(out0.q_log, 4); w_d = sig(out0.w_log, 3); p_d = sig(out0.p_log, 4);
    hw_d = jw * p_d;
    hb_d = j_body_eff * w_d' + wt * hw_d';
    h_inert = zeros(3, numel(t0));
    for k = 1:numel(t0)
        h_inert(:,k) = q2r(q_d(k,:)) * hb_d(:,k);
    end
    dh = vecnorm(h_inert - h_inert(:,1)) / norm(h_inert(:,1));
    fprintf('Max relative change in total inertial angular momentum |H(t)-H(0)|/|H(0)|: %.2e\n', max(dh));
    figure('Name', 'Sanity Check 1: Total Momentum Conservation');
    semilogy(t0, dh, 'LineWidth', 1.5); grid on;
    ylabel('|H_{tot}(t) - H_{tot}(0)| / |H_{tot}(0)|'); xlabel('Time (s)');
    title('Free-Drift Total Angular Momentum Conservation (Body + 4 Wheels)');
end

% =========================================================================
% RUN 1: IDEALIZED BLACK-BOX ACTUATION (ORIGINAL PAPER BASELINE)
% =========================================================================
fprintf('\n--- Run 1: Idealized Direct-Torque Actuation (Paper Baseline) ---\n');
set_param(g_ideal, 'Gain', '1'); set_param(g_rw, 'Gain', '0'); set_param(g_wdist, 'Gain', '0');
set_param(gt, 'Gain', num2str(target_on)); set_param(gd, 'Gain', num2str(dist_on));
out_ideal = sim(mdl, 'StopTime', num2str(tend));

t_id = out_ideal.q_log.time;
q_id = sig(out_ideal.q_log, 4);
w_id = sig(out_ideal.w_log, 3);
u_id = sig(out_ideal.u_log, 3);
tu_id = out_ideal.u_log.time;
qn_id = vecnorm(q_id(:,2:4), 2, 2);
k_id = find(qn_id > settle_tol*qn_id(1), 1, 'last');
tset_id = t_id(min(k_id + 1, numel(t_id)));
fprintf('Idealized Direct Torque: steady-state |q_v| = %.2e, 5%% settling time = %.1f s\n', qn_id(end), tset_id);

% =========================================================================
% RUN 2: REALISTIC TETRAHEDRAL REACTION WHEEL ACTUATION
% =========================================================================
fprintf('\n--- Run 2: Realistic Tetrahedral Reaction Wheel Actuation (%s) ---\n', rw_preset);
set_param(g_ideal, 'Gain', '0'); set_param(g_rw, 'Gain', '1'); set_param(g_wdist, 'Gain', num2str(noise_on));
set_param(gt, 'Gain', num2str(target_on)); set_param(gd, 'Gain', num2str(dist_on));
save_system(mdl); % Leave saved model in RW closed-loop mode
out_rw = sim(mdl, 'StopTime', num2str(tend));

t_rw = out_rw.q_log.time;
q_rw = sig(out_rw.q_log, 4);
w_rw = sig(out_rw.w_log, 3);
p_rw = sig(out_rw.p_log, 4);
mm_rw = sig(out_rw.mm_log, 4);
tm_rw = out_rw.mm_log.time;
hw_rw = jw * p_rw;

qn_rw = vecnorm(q_rw(:,2:4), 2, 2);
k_rw = find(qn_rw > settle_tol*qn_rw(1), 1, 'last');
tset_rw = t_rw(min(k_rw + 1, numel(t_rw)));
mmc = interp1(tm_rw, mm_rw, t_rw, 'previous', 'extrap');
energy_rw = trapz(t_rw, sum(abs(mmc .* p_rw), 2));
fprintf('Realistic RW Actuation:  steady-state |q_v| = %.2e, 5%% settling time = %.1f s\n', qn_rw(end), tset_rw);
fprintf('Actuator Mechanical Energy: %.3e J | Peak Wheel Momentum |h_w|/h_max: %.2f\n', ...
    energy_rw, max(abs(hw_rw(:)))/hmax);

% =========================================================================
% PLOTS: SANITY CHECK COMPARISON (IDEAL PAPER vs REALISTIC REACTION WHEELS)
% =========================================================================
figure('Name', 'Sanity Check 2: Idealized Paper vs Realistic Reaction Wheels', 'Position', [100 100 850 700]);
subplot(3,1,1);
plot(t_id, qn_id, 'b--', 'LineWidth', 1.8); hold on;
plot(t_rw, qn_rw, 'r-', 'LineWidth', 1.8); grid on;
yline(settle_tol*qn_id(1), 'k:', '5% Settling Threshold');
ylabel('Attitude Error ||q_v||_2');
title(sprintf('Sanity Check: Ideal Black-Box Torque (Settled: %.1fs) vs Realistic RWs (Settled: %.1fs)', tset_id, tset_rw));
legend('Idealized Black-Box (Paper, \tau_{max}=20 Nm)', sprintf('Tetrahedral RWs (\\tau_{w,max}=%.1f Nm)', taumax), 'Location', 'northeast');

subplot(3,1,2);
plot(t_id, vecnorm(w_id, 2, 2)*180/pi, 'b--', 'LineWidth', 1.8); hold on;
plot(t_rw, vecnorm(w_rw, 2, 2)*180/pi, 'r-', 'LineWidth', 1.8); grid on;
ylabel('Angular Rate ||\omega||_2 (deg/s)');
legend('Idealized Black-Box', 'Realistic RWs', 'Location', 'northeast');

subplot(3,1,3);
tau_rw_net = (-wt * mm_rw')'; % Net 3-axis torque delivered to body by wheels
stairs(tu_id, vecnorm(u_id, 2, 2), 'b--', 'LineWidth', 1.5); hold on;
stairs(tm_rw, vecnorm(tau_rw_net, 2, 2), 'r-', 'LineWidth', 1.5); grid on;
ylabel('Net Body Torque ||\tau||_2 (N m)'); xlabel('Time (s)');
legend('Ideal Commanded ||u||_2', 'Actual RW Delivered ||-W_T M_m||_2', 'Location', 'northeast');

% Detailed RW Telemetry Figure
figure('Name', 'Tetrahedral Reaction Wheel Telemetry', 'Position', [160 120 850 750]);
subplot(4,1,1);
plot(t_rw, q_rw(:,2:4), 'LineWidth', 1.4); grid on;
ylabel('q_v'); title('Closed-Loop Attitude Stabilization with Tetrahedral Reaction Wheels');
legend('q_1', 'q_2', 'q_3');

subplot(4,1,2);
plot(t_rw, w_rw*180/pi, 'LineWidth', 1.4); grid on;
ylabel('Body \omega (deg/s)'); legend('\omega_x', '\omega_y', '\omega_z');

subplot(4,1,3);
plot(t_rw, p_rw * 60 / (2*pi), 'LineWidth', 1.4); grid on;
ylabel('Wheel Speed (RPM)'); legend('RW 1', 'RW 2', 'RW 3', 'RW 4');

subplot(4,1,4);
stairs(tm_rw, mm_rw, 'LineWidth', 1.4); hold on;
yline([taumax -taumax], 'k--', 'Torque Limit'); grid on;
ylabel('Motor Torque M_m (N m)'); xlabel('Time (s)');
legend('RW 1', 'RW 2', 'RW 3', 'RW 4');

% ---------------- Helper Functions ----------------
function f = addrt(mdl, parent, name, d, r, xy)
b = add_block('sm_lib/Frames and Transforms/Rigid Transform', [mdl '/' name], 'Position', [xy xy + [50 40]]);
if ~isempty(r)
    sp(b, {'RotationMethod'}, {'RotationMatrix', 'Rotation Matrix'});
    sp(b, {'RotationMatrix'}, {mat2str(r, 17)});
end
if any(d)
    sp(b, {'TranslationMethod'}, {'Cartesian'});
    sp(b, {'TranslationCartesianOffset'}, {mat2str(d(:)', 17)});
end
add_line(mdl, parent, lc(b, 1));
f = rc(b, 1);
end

function f = addinertia(mdl, parent, name, mass, j, xy)
b = add_block('sm_lib/Body Elements/Inertia', [mdl '/' name], 'Position', [xy xy + [50 40]]);
sp(b, {'InertiaType'}, {'Custom', 'CustomInertia'});
sp(b, {'Mass'}, {num2str(mass, 17)});
sp(b, {'MomentsOfInertia'}, {mat2str(diag(j)', 17)});
sp(b, {'ProductsOfInertia'}, {mat2str([j(2,3) j(3,1) j(1,2)], 17)});
add_line(mdl, parent, one(b));
f = b;
end

function addbox_translucent(mdl, parent, name, ctr, dims, col, opacity, xy)
f = addrt(mdl, parent, [name ' Offset'], ctr, [], xy);
b = add_block('sm_lib/Body Elements/Brick Solid', [mdl '/' name], 'Position', [xy + [70 0] xy + [120 40]]);
sp(b, {'BrickDimensions'}, {mat2str(dims(:)', 17)});
sp(b, {'Density'}, {'1e-9'});
sp(b, {'GraphicDiffuseColor'}, {mat2str(col)}, true);
sp(b, {'GraphicOpacity'}, {num2str(opacity)}, true);
add_line(mdl, f, one(b));
end

function addcyl(mdl, parent, name, p1, p2, rad, col, xy)
a = (p2(:) - p1(:))/norm(p2 - p1);
r = [null(a') a];
if det(r) < 0, r(:,1) = -r(:,1); end
f = addrt(mdl, parent, [name ' Frame'], (p1(:) + p2(:))/2, r, xy);
b = add_block('sm_lib/Body Elements/Cylindrical Solid', [mdl '/' name], 'Position', [xy + [70 0] xy + [120 40]]);
sp(b, {'CylinderRadius'}, {num2str(rad, 17)});
sp(b, {'CylinderLength'}, {num2str(norm(p2 - p1), 17)});
sp(b, {'Density'}, {'1e-9'});
sp(b, {'GraphicDiffuseColor'}, {mat2str(col)}, true);
add_line(mdl, f, one(b));
end

function r = q2r(q)
s = q(1); v = q(2:4)';
r = (s^2 - v'*v)*eye(3) + 2*(v*v') + 2*s*[0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end

function ok = sp(b, names, vals, quiet)
if nargin < 4, quiet = false; end
for i = 1:numel(names)
    for j = 1:numel(vals)
        try
            set_param(b, names{i}, vals{j});
            ok = true;
            return
        catch
        end
    end
end
ok = false;
if ~quiet
    d = fieldnames(get_param(b, 'DialogParameters'));
    warning('%s: could not set %s. Its parameter names are:\n%s', getfullname(b), names{1}, strjoin(d', ', '));
end
end

function setscript(path, lines)
ch = sfroot().find('-isa', 'Stateflow.EMChart', 'Path', path);
ch.Script = strjoin(lines, newline);
end

function h = one(b)
p = get_param(b, 'PortHandles');
h = [p.LConn p.RConn];
h = h(1);
end

function h = lc(b, k), p = get_param(b, 'PortHandles'); h = p.LConn(k); end
function h = rc(b, k), p = get_param(b, 'PortHandles'); h = p.RConn(k); end
function h = ip(b, k), p = get_param(b, 'PortHandles'); h = p.Inport(k); end
function h = op(b, k), p = get_param(b, 'PortHandles'); h = p.Outport(k); end

function x = sig(s, n)
x = squeeze(s.signals.values);
if size(x, 2) ~= n, x = x.'; end
end