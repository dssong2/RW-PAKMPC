function rw = build_rw_array(mode)
%BUILD_RW_ARRAY  Tetrahedral 4-wheel Reaction Wheel Array block for any Simscape Multibody spacecraft.
%   rw = build_rw_array            builds rw_array_lib.slx with the block 'Reaction Wheel Array'
%   rw = build_rw_array('check')   also builds rw_array_check.slx and runs four checks
%   rw = build_rw_array('params')  same as rw_params(), kept for old calls
%   Parameter values are in rw_params.m.
%
%   Block: in  u (3x1 body torque command, N m), kn (null-space wheel torque, N m)
%          B   frame port, connect to the frame the wheels are mounted on (servicer body frame)
%          out p (wheel speed rel. body, rad/s), h = jw*p (N m s), tau = -W_T*M (N m),
%              sat (0 free, 1 torque limit, 2 momentum limit), pwr = M.*p (W),
%              M (applied motor torque per wheel after torque error and limits, N m)
%   The mask reads one struct, default name rw, so put rw = rw_params() in the workspace.
%   The rotors add rw.mass and rw.ilock (locked inertia about B) to the spacecraft: subtract them from the hull.
%
%   Dynamics: H = J_c w + J_w W_T p,  J_w(pdot_i + w_i'*wdot) = M_i - c_th p_i - t_c tanh(p_i/p_c)
%   Applied torque, once per controller sample: M = -W_T^+ u + kn n - (noise + bias + drift),
%   clipped to tau_max, and an accelerating wheel is held so |J_w p| cannot pass h_max within one sample.
if nargin < 1, mode = 'build'; end

here = fileparts(mfilename('fullpath'));
addpath(here);
rw = rw_params();                       % all wheel values live in rw_params.m
if strcmp(mode, 'params'), return, end

load_system('sm_lib'); load_system('nesl_utility'); load_system('simulink');
lib = 'rw_array_lib';
if bdIsLoaded('rw_array_check'), close_system('rw_array_check', 0); end
if bdIsLoaded(lib), close_system(lib, 0); end
new_system(lib, 'Library');
blk = [lib '/Reaction Wheel Array'];
add_block('simulink/Ports & Subsystems/Subsystem', blk, 'Position', [100 100 260 260]);
delete_line(blk, 'In1/1', 'Out1/1');
pos = @(x, y) [x y x+50 y+40];

set_param([blk '/In1'], 'Name', 'u', 'Position', pos(20, 40));
uin = [blk '/u'];
knin = add_block('simulink/Sources/In1', [blk '/kn'], 'Position', pos(20, 120));
bport = add_block('nesl_utility/Connection Port', [blk '/B'], 'Position', pos(20, 400));
sp(bport, {'Side'}, {'Left', 'left'}, true);
set_param([blk '/Out1'], 'Name', 'p', 'Position', pos(1300, 40));
outs = {'p', 'h', 'tau', 'sat', 'pwr', 'M'};
for k = 2:6
    add_block('simulink/Sinks/Out1', [blk '/' outs{k}], 'Position', pos(1300, 40 + 80*(k-1)));
end

% ---------------- wheels ----------------
mux = add_block('simulink/Signal Routing/Mux', [blk '/pMux'], 'Inputs', '4', 'Position', [900 400 910 900]);
dmx = add_block('simulink/Signal Routing/Demux', [blk '/MDemux'], 'Outputs', '4', 'Position', [660 400 670 900]);
for i = 1:4
    y = 400 + 130*(i-1);
    r = [null(rw.wt(:,i)') rw.wt(:,i)];
    if det(r) < 0, r(:,1) = -r(:,1); end
    rt = add_block('sm_lib/Frames and Transforms/Rigid Transform', sprintf('%s/Mount %d', blk, i), 'Position', pos(120, y));
    sp(rt, {'RotationMethod'}, {'RotationMatrix', 'Rotation Matrix'});
    sp(rt, {'RotationMatrix'}, {mat2str(r, 17)});
    sp(rt, {'TranslationMethod'}, {'Cartesian'});
    sp(rt, {'TranslationCartesianOffset'}, {sprintf('prm.d(:,%d)''', i)});

    rev = add_block('sm_lib/Joints/Revolute Joint', sprintf('%s/Wheel Joint %d', blk, i), 'Position', pos(240, y));
    sp(rev, {'TorqueActuationMode'}, {'InputTorque', 'Provided by Input'});
    sp(rev, {'MotionActuationMode'}, {'ComputedMotion', 'Automatically Computed'});
    sp(rev, {'SenseVelocity'}, {'on'});
    if sp(rev, {'DampingCoefficientUnits'}, {'N*m/(rad/s)'}, true)
        sp(rev, {'DampingCoefficient'}, {'prm.cth'});
    else
        sp(rev, {'DampingCoefficient'}, {'prm.cth*pi/180'});
    end

    whl = add_block('sm_lib/Body Elements/Inertia', sprintf('%s/Rotor %d', blk, i), 'Position', pos(360, y));
    sp(whl, {'InertiaType'}, {'Custom', 'CustomInertia'});
    sp(whl, {'Mass'}, {'prm.mw'});
    sp(whl, {'MomentsOfInertia'}, {'[prm.jt prm.jt prm.jw]'});

    % graphics only (near-massless), drawn where the rotor mass actually is
    cy = add_block('sm_lib/Body Elements/Cylindrical Solid', sprintf('%s/Rotor Graphic %d', blk, i), 'Position', pos(360, y+60));
    sp(cy, {'CylinderRadius'}, {'prm.rg'});
    sp(cy, {'CylinderLength'}, {'prm.lg'});
    sp(cy, {'Density'}, {'1e-9'});
    sp(cy, {'GraphicDiffuseColor'}, {'[0.9 0.45 0.1]'}, true);
    mo = add_block('sm_lib/Frames and Transforms/Rigid Transform', sprintf('%s/Marker Offset %d', blk, i), 'Position', pos(440, y+60));
    sp(mo, {'TranslationMethod'}, {'Cartesian'});
    sp(mo, {'TranslationCartesianOffset'}, {'[0.55*prm.rg 0 0]'});
    mb = add_block('sm_lib/Body Elements/Brick Solid', sprintf('%s/Marker %d', blk, i), 'Position', pos(520, y+60));
    sp(mb, {'BrickDimensions'}, {'[0.7*prm.rg 0.22*prm.rg 1.35*prm.lg]'});
    sp(mb, {'Density'}, {'1e-9'});
    sp(mb, {'GraphicDiffuseColor'}, {'[0.95 0.1 0.1]'}, true);

    pso = add_block('nesl_utility/PS-Simulink Converter', sprintf('%s/p%d', blk, i), 'Position', pos(800, y));
    sp(pso, {'Unit'}, {'rad/s'});
    spt = add_block('nesl_utility/Simulink-PS Converter', sprintf('%s/M%d', blk, i), 'Position', pos(720, y+60));
    sp(spt, {'Unit'}, {'N*m'});

    pr = get_param(rev, 'PortHandles');
    if numel(pr.RConn) < 2 || numel(pr.LConn) < 2
        error('Wheel Joint %d has %d left and %d right ports, expected 2 and 2.', i, numel(pr.LConn), numel(pr.RConn));
    end
    ln(blk, one(bport), lc(rt, 1));
    ln(blk, rc(rt, 1), pr.LConn(1));
    ln(blk, pr.RConn(1), one(whl));
    ln(blk, pr.RConn(1), one(cy));
    ln(blk, pr.RConn(1), lc(mo, 1));
    ln(blk, rc(mo, 1), one(mb));
    ln(blk, pr.RConn(2), lc(pso, 1));
    ln(blk, rc(spt, 1), pr.LConn(2));
    ln(blk, op(pso, 1), ip(mux, i));
    ln(blk, op(dmx, i), ip(spt, 1));
end

% ---------------- torque error: noise + bias + random-walk drift ----------------
nz = add_block('simulink/Sources/Random Number', [blk '/Noise'], 'Mean', '[0 0 0 0]', ...
    'Variance', 'max(prm.dist*prm.sigma^2, 1e-30)*ones(1,4)', 'Seed', 'prm.seed + (1:4)', ...
    'SampleTime', 'prm.ts', 'Position', pos(120, 200));
bs = add_block('simulink/Sources/Constant', [blk '/Bias'], 'Value', 'prm.dist*prm.b', 'Position', pos(120, 260));
dn = add_block('simulink/Sources/Random Number', [blk '/Drift Increment'], 'Mean', '[0 0 0 0]', ...
    'Variance', 'max(prm.dist*prm.sigd^2*prm.ts, 1e-30)*ones(1,4)', 'Seed', 'prm.seed + (11:14)', ...
    'SampleTime', 'prm.ts', 'Position', pos(120, 320));
di = add_block('simulink/Discrete/Discrete-Time Integrator', [blk '/Drift'], 'SampleTime', 'prm.ts', 'Position', pos(200, 320));
sp(di, {'gainval'}, {'1/prm.ts'});
ds = add_block('simulink/Math Operations/Sum', [blk '/Torque Error'], 'Inputs', '+++', 'Position', pos(280, 260));

% ---------------- allocation and limits, once per controller sample ----------------
zu = add_block('simulink/Discrete/Zero-Order Hold', [blk '/ZOHu'], 'SampleTime', 'prm.ts', 'Position', pos(120, 40));
zk = add_block('simulink/Discrete/Zero-Order Hold', [blk '/ZOHkn'], 'SampleTime', 'prm.ts', 'Position', pos(120, 120));
zp = add_block('simulink/Discrete/Zero-Order Hold', [blk '/ZOHp'], 'SampleTime', 'prm.ts', 'Position', pos(960, 300));
cwp = add_block('simulink/Sources/Constant', [blk '/W_T pinv'], 'Value', 'prm.wp', 'Position', pos(280, 40));
cn = add_block('simulink/Sources/Constant', [blk '/Null Vector'], 'Value', 'prm.n', 'Position', pos(280, 120));
cl = add_block('simulink/Sources/Constant', [blk '/Limit Values'], 'Value', '[prm.jw prm.taumax prm.hmax prm.ts]', 'Position', pos(280, 340));
lim = add_block('simulink/User-Defined Functions/MATLAB Function', [blk '/Limits'], 'Position', [400 40 540 360]);
setscript([blk '/Limits'], {
    'function [m, sat] = limits(u, kn, dw, p, wp, n, lim)'
    '% lim = [jw taumax hmax ts]'
    'm = -wp*u(:) + kn*n(:) - dw(:);'
    'p = p(:);'
    'sat = zeros(4, 1);'
    'for i = 1:4'
    '    if abs(m(i)) > lim(2)'
    '        m(i) = sign(m(i))*lim(2);'
    '        sat(i) = 1;'
    '    end'
    '    room = max(lim(3) - lim(1)*abs(p(i)), 0)/lim(4);'
    '    if sign(m(i)) == sign(p(i)) && abs(m(i)) > room'
    '        m(i) = sign(m(i))*room;'
    '        sat(i) = 2;'
    '    end'
    'end'}, {'m', 'sat'}, {[4 1], [4 1]});

% ---------------- Coulomb friction, continuous ----------------
cf = add_block('simulink/Sources/Constant', [blk '/Coulomb Values'], 'Value', '[prm.tc prm.pc]', 'Position', pos(960, 200));
fr = add_block('simulink/User-Defined Functions/MATLAB Function', [blk '/Coulomb'], 'Position', [1040 160 1140 240]);
setscript([blk '/Coulomb'], {
    'function tf = coulomb(p, c)'
    'tf = c(1)*tanh(p(:)/c(2));'}, {'tf'}, {[4 1]});
fs = add_block('simulink/Math Operations/Sum', [blk '/Net Motor Torque'], 'Inputs', '+-', 'Position', pos(580, 600));

% ---------------- outputs ----------------
gh = add_block('simulink/Math Operations/Gain', [blk '/Wheel Momentum'], 'Gain', 'prm.jw', 'Position', pos(1180, 120));
gt = add_block('simulink/Math Operations/Gain', [blk '/Body Torque'], 'Gain', '-prm.wt', 'Position', pos(1180, 200));
sp(gt, {'Multiplication'}, {'Matrix(K*u)'});
pw = add_block('simulink/Math Operations/Product', [blk '/Motor Power'], 'Inputs', '2', 'Position', pos(1180, 360));

ln(blk, op(uin, 1), ip(zu, 1));
ln(blk, op(knin, 1), ip(zk, 1));
ln(blk, op(nz, 1), ip(ds, 1));
ln(blk, op(bs, 1), ip(ds, 2));
ln(blk, op(dn, 1), ip(di, 1));
ln(blk, op(di, 1), ip(ds, 3));
ln(blk, op(mux, 1), ip(zp, 1));
ln(blk, op(zu, 1), ip(lim, 1));
ln(blk, op(zk, 1), ip(lim, 2));
ln(blk, op(ds, 1), ip(lim, 3));
ln(blk, op(zp, 1), ip(lim, 4));
ln(blk, op(cwp, 1), ip(lim, 5));
ln(blk, op(cn, 1), ip(lim, 6));
ln(blk, op(cl, 1), ip(lim, 7));
ln(blk, op(mux, 1), ip(fr, 1));
ln(blk, op(cf, 1), ip(fr, 2));
ln(blk, op(lim, 1), ip(fs, 1), 'M');
ln(blk, op(fr, 1), ip(fs, 2));
ln(blk, op(fs, 1), ip(dmx, 1));
ln(blk, op(mux, 1), ip([blk '/p'], 1), 'p');
ln(blk, op(mux, 1), ip(gh, 1));
ln(blk, op(gh, 1), ip([blk '/h'], 1));
ln(blk, op(lim, 1), ip(gt, 1));
ln(blk, op(gt, 1), ip([blk '/tau'], 1));
ln(blk, op(lim, 2), ip([blk '/sat'], 1));
ln(blk, op(lim, 1), ip(pw, 1));
ln(blk, op(mux, 1), ip(pw, 2));
ln(blk, op(pw, 1), ip([blk '/pwr'], 1));
ln(blk, op(lim, 1), ip([blk '/M'], 1));

msk = Simulink.Mask.create(blk);
msk.addParameter('Type', 'edit', 'Prompt', 'Parameter struct (from build_rw_array)', 'Name', 'prm', 'Value', 'rw');
msk.Description = ['Tetrahedral 4-wheel reaction wheel array. In: u body torque (3x1), kn null-space torque. ' ...
    'Out: p, h, tau, sat, pwr, M. Connect B to the mounting frame.'];
save_system(lib, fullfile(here, [lib '.slx']));
fprintf('Built %s.slx, block ''Reaction Wheel Array''. Wheel mass %.1f kg, locked inertia diag [%.4f %.4f %.4f] kg m^2.\n', ...
    lib, rw.mass, diag(rw.ilock));
if ~strcmp(mode, 'check'), return, end

% ================= check harness: free body + block =================
jc = [986.1 1 -1; 1 573.2 -1; -1 -1 1116];
jb = jc - rw.ilock;
if min(eig((jb + jb')/2)) <= 0, error('Hull inertia jc - rw.ilock is not positive definite.'); end
mdl = 'rw_array_check';
new_system(mdl);
set_param(mdl, 'Solver', 'ode45', 'RelTol', '1e-8', 'AbsTol', '1e-10', 'MaxStep', '0.01', 'StopTime', '10');
try, set_param(mdl, 'SimMechanicsOpenEditorOnUpdate', 'off'); catch, end
world = add_block('sm_lib/Frames and Transforms/World Frame', [mdl '/World'], 'Position', pos(30, 100));
mech = add_block('sm_lib/Utilities/Mechanism Configuration', [mdl '/Mech'], 'Position', pos(30, 180));
solver = add_block('nesl_utility/Solver Configuration', [mdl '/Solver'], 'Position', pos(30, 260));
if ~sp(mech, {'UniformGravity'}, {'None'}, true), sp(mech, {'GravityVector'}, {'[0 0 0]'}); end
six = add_block('sm_lib/Joints/6-DOF Joint', [mdl '/Body Joint'], 'Position', pos(150, 100));
sp(six, {'SphPositionTargetSpecify'}, {'on'});
sp(six, {'SphPositionTargetRotationMethod'}, {'Quaternion'});
sp(six, {'SphPositionTargetQuaternion'}, {'q0t'});
sp(six, {'SphVelocityTargetSpecify'}, {'on'});
sp(six, {'SphVelocityTargetValueUnits'}, {'rad/s'});
sp(six, {'SphVelocityTargetValue'}, {'w0t'});
sp(six, {'SphVelocityTargetInFollowerFrame'}, {'on'});
sp(six, {'SphSensePosition'}, {'on'});
sp(six, {'SphSenseVelocity'}, {'on'});
sp(six, {'SphSensingFrame'}, {'Follower', 'FollowerFrame', 'SensingFrameFollower'});
body = add_block('sm_lib/Body Elements/Inertia', [mdl '/Hull'], 'Position', pos(300, 20));
sp(body, {'InertiaType'}, {'Custom', 'CustomInertia'});
sp(body, {'Mass'}, {num2str(1155 - rw.mass, 17)});
sp(body, {'MomentsOfInertia'}, {mat2str(diag(jb)', 17)});
sp(body, {'ProductsOfInertia'}, {mat2str([jb(2,3) jb(3,1) jb(1,2)], 17)});
psq = add_block('nesl_utility/PS-Simulink Converter', [mdl '/Q'], 'Position', pos(300, 180));
psw = add_block('nesl_utility/PS-Simulink Converter', [mdl '/w'], 'Position', pos(300, 250));
sp(psq, {'Unit'}, {'1'}); sp(psw, {'Unit'}, {'rad/s'});
arr = add_block([lib '/Reaction Wheel Array'], [mdl '/Reaction Wheel Array'], 'Position', [500 60 660 260]);
cu = add_block('simulink/Sources/Constant', [mdl '/u'], 'Value', 'ucmd', 'Position', pos(330, 60));
ck = add_block('simulink/Sources/Constant', [mdl '/kn'], 'Value', 'kncmd', 'Position', pos(330, 120));

ph = get_param(six, 'PortHandles');
sq = [ph.RConn(2:end) ph.LConn(2:end)];
if numel(sq) ~= 2, error('Body Joint has %d sensing ports, expected 2.', numel(sq)); end
pa = get_param(arr, 'PortHandles');
ln(mdl, one(world), ph.LConn(1));
ln(mdl, one(world), one(mech));
ln(mdl, one(world), one(solver));
ln(mdl, ph.RConn(1), one(body));
ln(mdl, ph.RConn(1), pa.LConn(1));
ln(mdl, sq(1), lc(psq, 1));
ln(mdl, sq(2), lc(psw, 1));
ln(mdl, op(cu, 1), pa.Inport(1));
ln(mdl, op(ck, 1), pa.Inport(2));
logs = {'q_log', op(psq, 1); 'w_log', op(psw, 1); 'p_log', pa.Outport(1); 'tau_log', pa.Outport(3); 'sat_log', pa.Outport(4)};
for k = 1:size(logs, 1)
    lg = add_block('simulink/Sinks/To Workspace', [mdl '/' logs{k,1}], 'VariableName', logs{k,1}, ...
        'SaveFormat', 'Structure With Time', 'Position', pos(800, 40 + 70*k));
    ln(mdl, logs{k,2}, ip(lg, 1));
end
save_system(mdl, fullfile(here, [mdl '.slx']));

r0 = rw; r0.dist = 0;
q0 = [0.5 0.5 0.5 0.5]; w0 = [0.02 -0.01 -0.01];

s = runcheck(mdl, r0, q0, w0, [0;0;0], 0, 60);
hb = jc*s.w' + rw.jw*rw.wt*s.p';
hi = zeros(3, numel(s.t));
for k = 1:numel(s.t), hi(:,k) = q2r(s.q(k,:))*hb(:,k); end
fprintf('1 free drift, 60 s: max |H(t)-H(0)|/|H(0)| = %.2e (expect 1e-10 or less)\n', max(vecnorm(hi - hi(:,1)))/norm(hi(:,1)));

s = runcheck(mdl, r0, [1 0 0 0], [0 0 0], [0;0;0], 0.2, 20);
pe = s.p(end,:)';
fprintf('2 null space, kn = 0.2 N m, 20 s: max |w| = %.2e rad/s (expect ~0), p off n = %.2e (expect ~0), |p| = %.1f rad/s\n', ...
    max(vecnorm(s.w, 2, 2)), norm(pe - (rw.n'*pe)*rw.n)/norm(pe), norm(pe));

r3 = r0; r3.hmax = 1;
s = runcheck(mdl, r3, [1 0 0 0], [0 0 0], [0; 2; 0], 0, 10);
fprintf('3 momentum limit, h_max = 1 N m s: max |h|/h_max = %.4f (expect <= 1), wheels at limit at end: %d of 4, |tau(end)| = %.2e N m\n', ...
    max(abs(rw.jw*s.p(:)))/r3.hmax, sum(s.sat(end,:) == 2), norm(s.tau(end,:)));

r4 = r0; r4.cth = 0; r4.taumax = 1;     % keep u4 below the torque clip (u4 needs 0.27 N m on wheel 2)
u4 = [0.3; -0.2; 0.1];
s = runcheck(mdl, r4, [1 0 0 0], [0 0 0], u4, 0, 5);
wexp = (jc - rw.jw*(rw.wt*rw.wt'))\u4*5;
fprintf('4 torque delivery, 5 s, no friction: max |tau - u| = %.2e N m, |w(5) - w_expected|/|w_expected| = %.2e (expect ~1e-8)\n', ...
    max(vecnorm(s.tau - u4', 2, 2)), norm(s.w(end,:)' - wexp)/norm(wexp));
end

function s = runcheck(mdl, r, q0, w0, u, kn, tend)
in = Simulink.SimulationInput(mdl);
in = in.setVariable('rw', r);
in = in.setVariable('q0t', q0);
in = in.setVariable('w0t', w0);
in = in.setVariable('ucmd', u);
in = in.setVariable('kncmd', kn);
in = in.setModelParameter('StopTime', num2str(tend));
out = sim(in);
s.t = out.q_log.time;
s.q = sig(out.q_log, 4); s.w = sig(out.w_log, 3); s.p = sig(out.p_log, 4);
s.tau = sig(out.tau_log, 3); s.sat = sig(out.sat_log, 4);
end

function r = q2r(q)
a = q(1); v = q(2:4)';
r = (a^2 - v'*v)*eye(3) + 2*(v*v') + 2*a*[0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end

function x = sig(v, n)
x = squeeze(v.signals.values);
if size(x, 2) ~= n, x = x.'; end
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

function setscript(path, lines, outs, sizes)
ch = sfroot().find('-isa', 'Stateflow.EMChart', 'Path', path);
ch.Script = strjoin(lines, newline);
for k = 1:numel(outs)
    d = ch.find('-isa', 'Stateflow.Data', 'Name', outs{k});
    d.Props.Array.Size = mat2str(sizes{k});
end
end

function h = ln(sys, a, b, name)
try
    h = add_line(sys, a, b, 'autorouting', 'smart');
catch
    h = add_line(sys, a, b);
end
if nargin > 3, set_param(h, 'Name', name); end
end

function h = one(b)
p = get_param(b, 'PortHandles');
h = [p.LConn p.RConn];
h = h(1);
end

function h = lc(b, k), p = get_param(b, 'PortHandles'); h = p.LConn(k); end
function h = rc(b, k), p = get_param(b, 'PortHandles'); h = p.RConn(k); end
function h = ip(b, k)
if ischar(b), b = get_param(b, 'Handle'); end
p = get_param(b, 'PortHandles'); h = p.Inport(k);
end
function h = op(b, k)
if ischar(b), b = get_param(b, 'Handle'); end
p = get_param(b, 'PortHandles'); h = p.Outport(k);
end
