%BUILD_PAKMPC_MULTIBODY  Simscape Multibody model of the combined spacecraft in Sun et al.,
%   "Physics-based Online Adaptive Koopman MPC Attitude Control for Combined Spacecraft
%   with Dynamic Uncertainties," Sec. V-A and Fig. 5.
%   Builds pakmpc_multibody.slx next to this file, runs a free-drift check (no torques)
%   to confirm the assembled bodies add up to the paper's Jc, then runs the closed loop.
%
%   Standalone: uses a PD controller in place of PAKMPC, no other files needed.
%
%   [paper] = quoted from the paper
%   [guess] = not given in the paper, chosen here

noise_on = 1;               % Sec. V-C sensor noise
dist_on = 1;                % Sec. V-C external disturbance
target_on = 1;              % Eq. (59) target torque
drift_check = true;         % free-drift momentum check before the closed loop

v = ver;
if ~any(strcmp({v.Name}, 'Simscape')) || ~any(strcmp({v.Name}, 'Simscape Multibody'))
    error('Need Simscape and Simscape Multibody: Home > Add-Ons > Get Add-Ons.');
end

here = fileparts(mfilename('fullpath'));
if isempty(here) || contains(here, 'Editor_'), here = pwd; end

% ---------------- [paper] Sec. V-A and V-C ----------------
ms = 1080; mt = 75;                     % kg
js = diag([405 405 405]);               % kg m^2
jt = diag([36.8 37.5 36.8]);            % kg m^2
jc = [986.1 1 -1; 1 573.2 -1; -1 -1 1116];   % combined inertia, kg m^2
q0 = [0.5; 0.5; 0.5; 0.5];              % scalar first
w0 = [0.02; -0.01; -0.01];              % rad/s
umax = 20;                              % |u_i| <= 20 N m
ts = 0.1; tend = 200;                   % s
kpt = 0.02*jt; kdt = 0.05*jt;           % Eq. (59) target gains
noise_var = 1e-4;                       % N(0, 1e-4) on w and q_v
settle_tol = 0.05;                      % 5% error margin

% ---------------- geometry [guess] ----------------
ls = sqrt(6*js(1,1)/ms);                        % servicer cube side, J = m L^2/6 -> 1.5 m
lt = sqrt(6/mt*(trace(jt) - 2*diag(jt)'));      % target box matching Jt -> 1.73 x 1.70 x 1.73 m
gap = 0.2;                                      % docking ring length between the two faces, m
r = [0; ls/2 + gap + lt(2)/2; 0];               % target center in the servicer frame, m
c = mt*r/(ms + mt);                             % combined center of mass in the servicer frame
kp = 40; kd = 400;                              % PD gains in place of PAKMPC (zeta about 1)

% arms, panels and ring carry the rest of Jc, lumped at the combined CoM
pa = @(m, d) m*((d'*d)*eye(3) - d*(d'));
jst = js + pa(ms, -c) + jt + pa(mt, r - c);
jres = jc - jst;
e = sort(eig((jres + jres')/2));
fprintf('appendage inertia principal moments: [%.2f %.2f %.2f] kg m^2\n', e);
if e(1) <= 0 || e(3) > e(1) + e(2) + 1e-9
    error('Residual inertia Jc - (servicer + target) is not physical. Change r or gap.');
end

% ---------------- model ----------------
mdl = 'pakmpc_multibody';
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
    error(['Combined Body Joint ports: %d left, %d right (expected B, t and F, Q, w). ' ...
        'Its parameters starting with Sph are:\n%s\nPaste this to Claude.'], numel(ph.LConn), ...
        numel(ph.RConn), strjoin(f(strncmp(f, 'Sph', 3))', ', '));
end
add_line(mdl, one(world), ph.LConn(1));
add_line(mdl, one(world), one(mech));
add_line(mdl, one(world), one(solver));
cm = ph.RConn(1);           % combined center-of-mass frame

% ---------------- bodies ----------------
res = addinertia(mdl, cm, 'Appendages (arms, panels, ring)', 1e-6, jres, [280 40]);
sv = addrt(mdl, cm, 'Servicer Frame', -c, [], [280 120]);
addinertia(mdl, sv, 'Servicer', ms, js, [380 120]);
addbox(mdl, sv, 'Servicer Body', [0 0 0], [ls ls ls], [0.75 0.75 0.8], [460 120]);
tg = addrt(mdl, sv, 'Target Frame', r, [], [380 200]);
addinertia(mdl, tg, 'Target', mt, jt, [480 200]);
addbox(mdl, tg, 'Target Body', [0 0 0], lt, [0.85 0.7 0.3], [560 200]);

% docking ring between the faces
addcyl(mdl, sv, 'Docking Ring', [0 ls/2 0], [0 ls/2 + gap 0], 0.35, [0.3 0.3 0.3], [460 280]);

% two three-link arms from the servicer's +y face to the target's sides
k = 0;
for s = [-1 1]
    pts = [s*0.3*ls,       ls/2,        0;
           s*(lt(1)/2+0.35), ls/2 + 0.25, 0;
           s*(lt(1)/2+0.45), r(2) - 0.2,  0;
           s*lt(1)/2,        r(2),        0];
    for j = 1:3
        k = k + 1;
        addcyl(mdl, sv, sprintf('Arm %d Link %d', (s+3)/2, j), pts(j,:), pts(j+1,:), 0.05, ...
            [0.2 0.4 0.8], [460 + 120*j, 300 + 80*k]);
    end
end

% solar panels on the target's +-x sides
for s = [-1 1]
    i = (s+3)/2;
    addcyl(mdl, tg, sprintf('Panel Boom %d', i), [s*lt(1)/2 0 0], [s*(lt(1)/2 + 0.25) 0 0], 0.03, ...
        [0.5 0.5 0.5], [700, 120 + 80*i]);
    addbox(mdl, tg, sprintf('Solar Panel %d', i), [s*(lt(1)/2 + 0.25 + 0.9) 0 0], [1.8 0.9 0.03], ...
        [0.1 0.15 0.45], [820, 120 + 80*i]);
end

% ---------------- sensing ----------------
psq = add_block('nesl_utility/PS-Simulink Converter', [mdl '/Q'], 'Position', pos(280, 520));
psw = add_block('nesl_utility/PS-Simulink Converter', [mdl '/w'], 'Position', pos(280, 590));
sp(psq, {'Unit'}, {'1'});
sp(psw, {'Unit'}, {'rad/s'});
add_line(mdl, ph.RConn(2), lc(psq, 1));
add_line(mdl, ph.RConn(3), lc(psw, 1));

zq = add_block('simulink/Discrete/Zero-Order Hold', [mdl '/ZOHq'], 'SampleTime', num2str(ts), 'Position', pos(400, 520));
zw = add_block('simulink/Discrete/Zero-Order Hold', [mdl '/ZOHw'], 'SampleTime', num2str(ts), 'Position', pos(400, 590));
nz = add_block('simulink/Sources/Random Number', [mdl '/Noise'], 'Mean', m(zeros(1,6)), ...
    'Variance', m(noise_var*ones(1,6)), 'Seed', m(1:6), 'SampleTime', num2str(ts), 'Position', pos(400, 660));
gn = add_block('simulink/Math Operations/Gain', [mdl '/Noise On'], 'Gain', num2str(noise_on), 'Position', pos(480, 660));
sen = add_block('simulink/User-Defined Functions/MATLAB Function', [mdl '/Sensor'], 'Position', [560 520 680 680]);
setscript([mdl '/Sensor'], {
    'function [wm, qm] = sensor(q, w, n)'
    '% Sec. V-C: N(0,1e-4) noise on w and q_v, quaternion renormalized'
    'q = q(:); w = w(:); n = n(:);'
    'wm = w + n(1:3);'
    'qm = [q(1); q(2:4) + n(4:6)];'
    'qm = qm/norm(qm);'});

% ---------------- controller ----------------
ctl = add_block('simulink/User-Defined Functions/MATLAB Function', [mdl '/Controller'], 'Position', [740 520 880 640]);
setscript([mdl '/Controller'], {
    'function tau = pd(wm, qm)'
    '% stand-in for PAKMPC: drive q to [1 0 0 0] and w to 0'
    sprintf('kp = %.17g; kd = %.17g;', kp, kd)
    's = 1; if qm(1) < 0, s = -1; end'
    'tau = -kp*s*qm(2:4) - kd*wm;'});
sat = add_block('simulink/Discontinuities/Saturation', [mdl '/Torque Limit'], 'UpperLimit', num2str(umax), ...
    'LowerLimit', num2str(-umax), 'Position', pos(920, 540));
gc = add_block('simulink/Math Operations/Gain', [mdl '/Control On'], 'Gain', '1', 'Position', pos(990, 540));

% ---------------- target torque, Eq. (59), and disturbance, Sec. V-C ----------------
tgt = add_block('simulink/User-Defined Functions/MATLAB Function', [mdl '/Target Torque'], 'Position', [560 740 680 800]);
setscript([mdl '/Target Torque'], {
    'function ut = target(q, w)'
    '% u_t = -Kpt q_vt - Kdt w, target holds its initial attitude (README: hold_initial)'
    'q = q(:); w = w(:);'
    ['a = ' m([q0(1); -q0(2:4)]) ';']
    ['kpt = ' m(kpt) '; kdt = ' m(kdt) ';']
    'qv = a(1)*q(2:4) + q(1)*a(2:4) + cross(a(2:4), q(2:4));'
    'ut = -kpt*qv - kdt*w;'});
gt = add_block('simulink/Math Operations/Gain', [mdl '/Target On'], 'Gain', num2str(target_on), 'Position', pos(720, 750));

clk = add_block('simulink/Sources/Clock', [mdl '/Clock'], 'Position', pos(480, 840));
dst = add_block('simulink/User-Defined Functions/MATLAB Function', [mdl '/Disturbance'], 'Position', [560 830 680 880]);
setscript([mdl '/Disturbance'], {
    'function td = dist(t)'
    'td = [0.1*sin(0.075*t); -0.2*sin(0.1*t); 0.15*sin(0.1*t)];'});
gd = add_block('simulink/Math Operations/Gain', [mdl '/Disturbance On'], 'Gain', num2str(dist_on), 'Position', pos(720, 840));

tsum = add_block('simulink/Math Operations/Sum', [mdl '/Torque Sum'], 'Inputs', '+++', 'Position', pos(1060, 640));
spt = add_block('nesl_utility/Simulink-PS Converter', [mdl '/Torque'], 'Position', pos(1140, 640));
sp(spt, {'Unit'}, {'N*m'});

add_line(mdl, op(psq, 1), ip(zq, 1));
add_line(mdl, op(psw, 1), ip(zw, 1));
add_line(mdl, op(nz, 1), ip(gn, 1));
add_line(mdl, op(zq, 1), ip(sen, 1));
add_line(mdl, op(zw, 1), ip(sen, 2));
add_line(mdl, op(gn, 1), ip(sen, 3));
add_line(mdl, op(sen, 1), ip(ctl, 1));
add_line(mdl, op(sen, 2), ip(ctl, 2));
add_line(mdl, op(ctl, 1), ip(sat, 1));
add_line(mdl, op(sat, 1), ip(gc, 1));
add_line(mdl, op(psq, 1), ip(tgt, 1));
add_line(mdl, op(psw, 1), ip(tgt, 2));
add_line(mdl, op(tgt, 1), ip(gt, 1));
add_line(mdl, op(clk, 1), ip(dst, 1));
add_line(mdl, op(dst, 1), ip(gd, 1));
add_line(mdl, op(gc, 1), ip(tsum, 1));
add_line(mdl, op(gt, 1), ip(tsum, 2));
add_line(mdl, op(gd, 1), ip(tsum, 3));
add_line(mdl, op(tsum, 1), ip(spt, 1));
add_line(mdl, rc(spt, 1), ph.LConn(2));

logs = {'q_log', psq, 1; 'w_log', psw, 1; 'tau_log', gc, 1; 'ut_log', gt, 1};
for i = 1:size(logs, 1)
    lg = add_block('simulink/Sinks/To Workspace', [mdl '/' logs{i,1}], 'VariableName', logs{i,1}, ...
        'SaveFormat', 'Structure With Time', 'Position', pos(1250, 40 + 70*i));
    add_line(mdl, op(logs{i,2}, logs{i,3}), ip(lg, 1));
end

save_system(mdl, fullfile(here, [mdl '.slx']));

% ---------------- free-drift check: no torques, H = R(q) Jc w must stay constant ----------------
if drift_check
    set_param(gc, 'Gain', '0'); set_param(gt, 'Gain', '0'); set_param(gd, 'Gain', '0');
    out = sim(mdl, 'StopTime', '100');
    t = out.q_log.time; q = sig(out.q_log, 4); w = sig(out.w_log, 3);
    h = zeros(3, numel(t));
    for k = 1:numel(t)
        h(:,k) = q2r(q(k,:))*jc*w(k,:)';
    end
    dh = vecnorm(h - h(:,1))/norm(h(:,1));
    fprintf('initial q from Simscape: [%.4f %.4f %.4f %.4f] (paper [0.5 0.5 0.5 0.5])\n', q(1,:));
    fprintf('initial w from Simscape: [%.4f %.4f %.4f] rad/s (paper [0.02 -0.01 -0.01])\n', w(1,:));
    fprintf('free drift, max |H(t) - H(0)|/|H(0)| using the paper Jc: %.2e\n', max(dh));
    figure('Name', 'Free-drift momentum check'); semilogy(t, dh); grid on;
    ylabel('|H(t) - H(0)| / |H(0)|'); xlabel('Time (s)');
    set_param(gc, 'Gain', '1'); set_param(gt, 'Gain', num2str(target_on)); set_param(gd, 'Gain', num2str(dist_on));
    save_system(mdl);
end

% ---------------- closed loop ----------------
out = sim(mdl);
t = out.q_log.time; q = sig(out.q_log, 4); w = sig(out.w_log, 3);
tau = sig(out.tau_log, 3);
qn = vecnorm(q(:,2:4), 2, 2);
k = find(qn > settle_tol*qn(1), 1, 'last');
fprintf('PD: steady-state |q_v| = %.2e, settled (5%%) at %.1f s\n', qn(end), t(min(k + 1, numel(t))));

figure('Name', 'Closed loop (PD)');
subplot(3,1,1); plot(t, q(:,2:4)); grid on; ylabel('q_v'); legend('q_1', 'q_2', 'q_3');
subplot(3,1,2); plot(t, w); grid on; ylabel('\omega (rad/s)'); legend('\omega_x', '\omega_y', '\omega_z');
subplot(3,1,3); stairs(out.tau_log.time, tau); grid on; ylabel('\tau (N m)'); xlabel('Time (s)');

% ---------------- helpers ----------------
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

function addbox(mdl, parent, name, ctr, dims, col, xy)
f = addrt(mdl, parent, [name ' Offset'], ctr, [], xy);
b = add_block('sm_lib/Body Elements/Brick Solid', [mdl '/' name], 'Position', [xy + [70 0] xy + [120 40]]);
sp(b, {'BrickDimensions'}, {mat2str(dims(:)', 17)});
sp(b, {'Density'}, {'1e-9'});
sp(b, {'GraphicDiffuseColor'}, {mat2str(col)}, true);
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
