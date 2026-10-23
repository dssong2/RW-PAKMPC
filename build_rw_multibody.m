%BUILD_RW_MULTIBODY  Rigid body with the tetrahedral 4-wheel reaction wheel model of the proposal,
%   "Effect of Reaction Wheel Actuator Dynamics on Adaptive Koopman MPC for Space Debris Capture," Sec. IV.
%   Builds rw_sim.slx next to this file, runs a PD law as a stand-in for PAKMPC's u*_0|k,
%   checks that total angular momentum is conserved, and prints the Sec. V metrics.

v = ver;
if ~any(strcmp({v.Name}, 'Simscape')) || ~any(strcmp({v.Name}, 'Simscape Multibody'))
    error('Need Simscape and Simscape Multibody: Home > Add-Ons > Get Add-Ons.');
end

% ---------------- spacecraft (arbitrary) ----------------
ib = diag([12 10 8]);       % kg m^2
mb = 100;                   % kg
w0 = [0.05; -0.03; 0.04];   % rad/s

% ---------------- reaction wheels, Sec. IV ----------------
beta = acos(1/sqrt(3));                          % Eq. (1), 54.74 deg
wt = [sin(beta) -sin(beta) 0          0;
      cos(beta)  cos(beta) cos(beta)  cos(beta);
      0          0         sin(beta) -sin(beta)];  % W_T, column i is the spin axis of wheel i
jw = 0.01;              % J_w = J in Eq. (3), rotor inertia about the shaft, kg m^2 [assumed]
jt = jw/2;              % rotor transverse inertia (thin disk), not in the proposal [assumed]
mw = 0.5;               % rotor mass, kg [assumed]
cth = 1e-5;             % c_theta in Eq. (3), bearing viscous drag, N m s/rad [assumed]
taumax = 0.1;           % |tau_w,i| <= tau_max, N m [assumed]
hmax = 4;               % |h_w,i| <= h_max, N m s [assumed]
sigma = 1e-4;           % Eq. (5) noise N(0, sigma^2 I4), N m [assumed]
b = 2e-4*[1 -1 0.5 -0.5];   % Eq. (5) bias b(mu,T), constant per wheel, N m [assumed]
sigd = 1e-5;            % Eq. (5) drift d(t), random walk, N m per sqrt(s) [assumed]

% ---------------- controller (stand-in for PAKMPC) ----------------
kp = 1; kd = 4;
ts = 0.1; tend = 120;
wp = pinv(wt);

% ---------------- model ----------------
mdl = 'rw_sim';
here = fileparts(mfilename('fullpath'));
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

six = add_block('sm_lib/Joints/6-DOF Joint', [mdl '/Body Joint'], 'Position', pos(150, 300));
f = fieldnames(get_param(six, 'DialogParameters'));
sp(six, {'SphVelocityTargetSpecify'}, {'on'});
sp(six, {'SphVelocityTargetValueUnits'}, {'rad/s'});
sp(six, {'SphVelocityTargetValue'}, {m(w0')});
sp(six, {'SphVelocityTargetInFollowerFrame'}, {'on'});
sp(six, {'SphSensePosition'}, {'on'});
sp(six, {'SphSenseVelocity'}, {'on'});
sp(six, {'SphSensingFrame'}, {'Follower', 'FollowerFrame', 'SensingFrameFollower'});

body = add_block('sm_lib/Body Elements/Inertia', [mdl '/Body'], 'Position', pos(280, 200));
sp(body, {'InertiaType'}, {'Custom', 'CustomInertia'});
sp(body, {'Mass'}, {num2str(mb)});
sp(body, {'MomentsOfInertia'}, {m(diag(ib)')});
sp(body, {'ProductsOfInertia'}, {m([ib(2,3) ib(3,1) ib(1,2)])});

% graphics only: near-massless box sized to match ib
dims = sqrt(6/mb*(trace(ib) - 2*diag(ib)'));
box = add_block('sm_lib/Body Elements/Brick Solid', [mdl '/Body Graphic'], 'Position', pos(280, 120));
sp(box, {'BrickDimensions'}, {m(dims)});
sp(box, {'Density'}, {'1e-9'});
sp(box, {'GraphicDiffuseColor'}, {'[0.6 0.65 0.75]'}, true);

psq = add_block('nesl_utility/PS-Simulink Converter', [mdl '/Q'], 'Position', pos(280, 330));
psw = add_block('nesl_utility/PS-Simulink Converter', [mdl '/w'], 'Position', pos(280, 400));
sp(psq, {'Unit'}, {'1'});
sp(psw, {'Unit'}, {'rad/s'});

mux = add_block('simulink/Signal Routing/Mux', [mdl '/pMux'], 'Inputs', '4', 'Position', [820 520 830 720]);
dmx = add_block('simulink/Signal Routing/Demux', [mdl '/MmDemux'], 'Outputs', '4', 'Position', [1250 520 1260 720]);

ph = get_param(six, 'PortHandles');
sq = [ph.RConn(2:end) ph.LConn(2:end)];
if numel(sq) ~= 2
    error(['Body Joint has %d sensing ports, expected 2 (Q, w). Its parameters with "sens" in the name are:\n%s\n' ...
        'Paste this to Claude.'], numel(sq), strjoin(f(~cellfun(@isempty, regexpi(f, 'sens')))', ', '));
end
add_line(mdl, one(world), ph.LConn(1));
add_line(mdl, one(world), one(mech));
add_line(mdl, one(world), one(solver));
add_line(mdl, ph.RConn(1), one(body));
add_line(mdl, ph.RConn(1), one(box));
add_line(mdl, sq(1), lc(psq, 1));
add_line(mdl, sq(2), lc(psw, 1));
bodyf = ph.RConn(1);

% ---------------- wheels: J_w pdot + c_theta p = M_m, spin axis = column i of W_T ----------------
for i = 1:4
    y = 500 + 110*(i-1);
    r = [null(wt(:,i)') wt(:,i)];
    if det(r) < 0, r(:,1) = -r(:,1); end

    rt = add_block('sm_lib/Frames and Transforms/Rigid Transform', sprintf('%s/RT%d', mdl, i), 'Position', pos(280, y));
    sp(rt, {'RotationMethod'}, {'RotationMatrix', 'Rotation Matrix'});
    sp(rt, {'RotationMatrix'}, {m(r)});

    rev = add_block('sm_lib/Joints/Revolute Joint', sprintf('%s/Wheel Joint %d', mdl, i), 'Position', pos(400, y));
    sp(rev, {'TorqueActuationMode'}, {'InputTorque', 'Provided by Input'});
    sp(rev, {'MotionActuationMode'}, {'ComputedMotion', 'Automatically Computed'});
    sp(rev, {'SenseVelocity'}, {'on'});
    if sp(rev, {'DampingCoefficientUnits'}, {'N*m/(rad/s)'}, true)
        sp(rev, {'DampingCoefficient'}, {num2str(cth, 17)});
    else
        sp(rev, {'DampingCoefficient'}, {num2str(cth*pi/180, 17)});    % dialog default is N*m/(deg/s)
    end

    whl = add_block('sm_lib/Body Elements/Inertia', sprintf('%s/Wheel %d', mdl, i), 'Position', pos(540, y));
    sp(whl, {'InertiaType'}, {'Custom', 'CustomInertia'});
    sp(whl, {'Mass'}, {num2str(mw)});
    sp(whl, {'MomentsOfInertia'}, {m([jt jt jw])});

    pso = add_block('nesl_utility/PS-Simulink Converter', sprintf('%s/p%d', mdl, i), 'Position', pos(680, y));
    sp(pso, {'Unit'}, {'rad/s'});
    spt = add_block('nesl_utility/Simulink-PS Converter', sprintf('%s/Mm%d', mdl, i), 'Position', pos(1300, y));
    sp(spt, {'Unit'}, {'N*m'});

    pr = get_param(rev, 'PortHandles');
    if numel(pr.RConn) < 2 || numel(pr.LConn) < 2
        fw = fieldnames(get_param(rev, 'DialogParameters'));
        error('Wheel Joint %d ports: %d left, %d right (expected B, t and F, w). Its parameters are:\n%s\nPaste this to Claude.', ...
            i, numel(pr.LConn), numel(pr.RConn), strjoin(fw', ', '));
    end
    add_line(mdl, bodyf, lc(rt, 1));
    add_line(mdl, rc(rt, 1), pr.LConn(1));
    add_line(mdl, pr.RConn(1), one(whl));

    % graphics only: near-massless disk pushed out along the spin axis so it shows outside the box
    gt = add_block('sm_lib/Frames and Transforms/Rigid Transform', sprintf('%s/Wheel Graphic Offset %d', mdl, i), 'Position', pos(540, y-50));
    sp(gt, {'TranslationMethod'}, {'StandardAxis', 'Standard Axis'});
    sp(gt, {'TranslationStandardAxis'}, {'+Z', 'PosZ'});
    sp(gt, {'TranslationStandardOffset'}, {num2str(0.6*max(dims))});
    cyl = add_block('sm_lib/Body Elements/Cylindrical Solid', sprintf('%s/Wheel Graphic %d', mdl, i), 'Position', pos(620, y-50));
    sp(cyl, {'CylinderRadius'}, {num2str(0.2*max(dims))});
    sp(cyl, {'CylinderLength'}, {num2str(0.05*max(dims))});
    sp(cyl, {'Density'}, {'1e-9'});
    sp(cyl, {'GraphicDiffuseColor'}, {'[0.9 0.45 0.1]'}, true);
    add_line(mdl, pr.RConn(1), lc(gt, 1));
    add_line(mdl, rc(gt, 1), one(cyl));

    add_line(mdl, pr.RConn(2), lc(pso, 1));
    add_line(mdl, op(pso, 1), ip(mux, i));
    add_line(mdl, op(dmx, i), ip(spt, 1));
    add_line(mdl, rc(spt, 1), pr.LConn(2));
end

% ---------------- controller ----------------
zq = add_block('simulink/Discrete/Zero-Order Hold', [mdl '/ZOHq'], 'SampleTime', num2str(ts), 'Position', pos(400, 330));
zw = add_block('simulink/Discrete/Zero-Order Hold', [mdl '/ZOHw'], 'SampleTime', num2str(ts), 'Position', pos(400, 400));
zp = add_block('simulink/Discrete/Zero-Order Hold', [mdl '/ZOHp'], 'SampleTime', num2str(ts), 'Position', pos(870, 600));

ctrl = add_block('simulink/User-Defined Functions/MATLAB Function', [mdl '/Controller'], 'Position', [500 330 620 430]);
setscript([mdl '/Controller'], {
    'function u = ctrl(q, w)'
    '% stand-in for PAKMPC u*_0|k'
    'q = q(:); w = w(:);'
    sprintf('kp = %.17g; kd = %.17g;', kp, kd)
    's = 1; if q(1) < 0, s = -1; end'
    'u = -kp*s*q(2:4) - kd*w;'});

% ---------------- Eq. (5) disturbances on each wheel: N(0, sigma^2 I4) + b + d(t) ----------------
nz = add_block('simulink/Sources/Random Number', [mdl '/Noise'], 'Mean', m(zeros(1,4)), ...
    'Variance', m(sigma^2*ones(1,4)), 'Seed', m(1:4), 'SampleTime', num2str(ts), 'Position', pos(800, 200));
bs = add_block('simulink/Sources/Constant', [mdl '/Bias'], 'Value', m(b), 'Position', pos(800, 260));
dn = add_block('simulink/Sources/Random Number', [mdl '/Drift Increment'], 'Mean', m(zeros(1,4)), ...
    'Variance', m(sigd^2*ts*ones(1,4)), 'Seed', m(11:14), 'SampleTime', num2str(ts), 'Position', pos(720, 320));
di = add_block('simulink/Discrete/Discrete-Time Integrator', [mdl '/Drift'], 'SampleTime', num2str(ts), ...
    'Position', pos(800, 320));
sp(di, {'gainval'}, {num2str(1/ts, 17)});    % d_k = sum of increments (random walk)
ds = add_block('simulink/Math Operations/Sum', [mdl '/Wheel Disturbance'], 'Inputs', '+++', 'Position', pos(900, 260));

% ---------------- allocation and limits ----------------
alloc = add_block('simulink/User-Defined Functions/MATLAB Function', [mdl '/Allocation'], 'Position', [1000 330 1140 430]);
setscript([mdl '/Allocation'], {
    'function mm = alloc(u, dw, p)'
    '% pseudo-inverse of W_T; body torque is -W_T*mm = u + W_T*dw, Eqs. (5)-(6)'
    'u = u(:); dw = dw(:); p = p(:);'
    ['wp = ' m(wp) ';']
    sprintf('jw = %.17g; taumax = %.17g; hmax = %.17g;', jw, taumax, hmax)
    'mm = -wp*u - dw;'
    'mm = min(max(mm, -taumax), taumax);'
    'mm(abs(jw*p) >= hmax & sign(mm) == sign(p)) = 0;'});

add_line(mdl, op(psq, 1), ip(zq, 1));
add_line(mdl, op(psw, 1), ip(zw, 1));
add_line(mdl, op(mux, 1), ip(zp, 1));
add_line(mdl, op(zq, 1), ip(ctrl, 1));
add_line(mdl, op(zw, 1), ip(ctrl, 2));
add_line(mdl, op(nz, 1), ip(ds, 1));
add_line(mdl, op(bs, 1), ip(ds, 2));
add_line(mdl, op(dn, 1), ip(di, 1));
add_line(mdl, op(di, 1), ip(ds, 3));
add_line(mdl, op(ctrl, 1), ip(alloc, 1));
add_line(mdl, op(ds, 1), ip(alloc, 2));
add_line(mdl, op(zp, 1), ip(alloc, 3));
add_line(mdl, op(alloc, 1), ip(dmx, 1));

logs = {'q_log', psq; 'w_log', psw; 'p_log', mux; 'mm_log', alloc; 'u_log', ctrl};
for i = 1:size(logs, 1)
    lg = add_block('simulink/Sinks/To Workspace', [mdl '/' logs{i,1}], 'VariableName', logs{i,1}, ...
        'SaveFormat', 'Structure With Time', 'Position', pos(1400, 40 + 70*i));
    add_line(mdl, op(logs{i,2}, 1), ip(lg, 1));
end

save_system(mdl, fullfile(here, [mdl '.slx']));

% ---------------- run ----------------
out = sim(mdl);
t = out.q_log.time;
q = sig(out.q_log, 4); w = sig(out.w_log, 3); p = sig(out.p_log, 4);
tm = out.mm_log.time; mm = sig(out.mm_log, 4);
hw = jw*p;

% ---------------- momentum conservation check ----------------
itot = ib;
for i = 1:4
    itot = itot + jt*(eye(3) - wt(:,i)*wt(:,i)') + jw*(wt(:,i)*wt(:,i)');
end
hb = itot*w' + wt*hw';
h = zeros(3, numel(t));
for k = 1:numel(t)
    s = q(k,1); qv = q(k,2:4)';
    r = (s^2 - qv'*qv)*eye(3) + 2*(qv*qv') + 2*s*[0 -qv(3) qv(2); qv(3) 0 -qv(1); -qv(2) qv(1) 0];
    h(:,k) = r*hb(:,k);
end
dh = vecnorm(h - h(:,1)) / norm(h(:,1));

% ---------------- Sec. V metrics ----------------
qn = vecnorm(q(:,2:4), 2, 2);
k = find(qn > 1e-3, 1, 'last');
tconv = t(min(k + 1, numel(t)));
mmc = interp1(tm, mm, t, 'previous', 'extrap');
energy = trapz(t, sum(abs(mmc.*p), 2));
fprintf('max relative change in inertial angular momentum: %.2e\n', max(dh));
fprintf('steady-state |q_v|: %.2e   time to |q_v| < 1e-3: %.1f s\n', qn(end), tconv);
fprintf('actuator energy (integral of sum |M_m,i p_i| dt): %.3e J\n', energy);
fprintf('peak |h_w,i| / h_max: %.2f\n', max(abs(hw(:))) / hmax);

% ---------------- plots ----------------
figure('Name', 'RW detumble');
subplot(4,1,1); plot(t, q(:,2:4)); grid on; ylabel('q_v'); legend('q_1', 'q_2', 'q_3');
subplot(4,1,2); plot(t, w*180/pi); grid on; ylabel('body \omega (deg/s)');
subplot(4,1,3); plot(t, p); grid on; ylabel('wheel p (rad/s)'); legend('1', '2', '3', '4');
subplot(4,1,4); stairs(tm, mm); grid on; ylabel('M_m (N m)'); xlabel('Time (s)');
figure('Name', 'Body rates');
plot(t, w*180/pi); grid on; ylabel('body \omega (deg/s)'); xlabel('Time (s)'); legend('\omega_x', '\omega_y', '\omega_z');
figure('Name', 'Wheel momentum');
plot(t, hw); hold on; yline([-hmax hmax], 'k--'); grid on;
ylabel('h_{w,i} = J_w p_i (N m s)'); xlabel('Time (s)'); legend('1', '2', '3', '4');
figure('Name', 'Momentum check'); semilogy(t, dh); grid on; ylabel('|H(t) - H(0)| / |H(0)|'); xlabel('Time (s)');

% ---------------- helpers ----------------
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
