function out = pakmpc(w0, q0, varargin)
%PAKMPC  Run PAKMPC (or a baseline) on the combined spacecraft from a given initial state.
%   out = pakmpc(w0, q0) runs PAKMPC from angular velocity w0 (rad/s, 3x1) and attitude q0
%   (scalar-first quaternion, normalized here) under the Sec. V-C scenario: nominal inertia
%   Jc0 = Js in the controller, true Jc in the ode45 plant, target torque, disturbance and noise.
%   pakmpc() uses the paper's initial state, w0 = [0.02 -0.01 -0.01], q0 = [0.5 0.5 0.5 0.5].
%
%   out = pakmpc(w0, q0, Name, Value) with options:
%     'controller'  'pakmpc' (default), 'nonadaptive', 'kmae' (KMAE-AKMPC) or 'inmpc' (INMPC)
%     'tend'        final time in s, default 200
%     'jc'          true inertia of the plant, default the paper's Jc
%     'dist'        include the disturbance tau_d, default true
%     'target'      include the target torque of Eq. (59), default true
%     'noise'       add measurement noise, default true
%     'seed'        noise seed, default 1
%     'params'      a pakmpc_params struct to use instead of the defaults
%     'plot'        plot q_v, w, tau and (PAKMPC) mu and ||Delta B||, default true
%
%   out fields: t (1 x nt), x = [w; q] (7 x nt), tau (3 x nt-1), mu, db (||Delta B||),
%   steptime (controller time per sample), settle (5% settling time, s), ss_qv and ss_w
%   (mean norms over the last 10 s), effort (integral of ||tau||^2).
%
%   Examples
%     out = pakmpc;                                       % the paper's case
%     out = pakmpc([0; 0.01; 0], [1; 0.3; -0.2; 0.1]);    % another initial state
%     out = pakmpc('controller', 'nonadaptive', 'noise', false, 'tend', 100);   % paper's state

if nargin >= 1 && (ischar(w0) || isstring(w0))     % options only: pakmpc('Name', Value, ...)
    if nargin >= 2, varargin = [{w0, q0} varargin]; else, varargin = {w0}; end
    w0 = []; q0 = [];
elseif nargin < 2
    q0 = [];
end
if nargin < 1, w0 = []; end

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'lib'));

ip = inputParser;
ip.addParameter('controller', 'pakmpc', @(s) any(strcmpi(s, {'pakmpc', 'nonadaptive', 'kmae', 'inmpc'})));
ip.addParameter('tend', [], @(v) isempty(v) || (isscalar(v) && v > 0));
ip.addParameter('jc', [], @(v) isempty(v) || isequal(size(v), [3 3]));
ip.addParameter('dist', true, @islogical);
ip.addParameter('target', true, @islogical);
ip.addParameter('noise', true, @islogical);
ip.addParameter('seed', 1, @isscalar);
ip.addParameter('params', [], @(v) isempty(v) || isstruct(v));
ip.addParameter('plot', true, @islogical);
ip.parse(varargin{:});
o = ip.Results;

p = o.params;
if isempty(p), p = pakmpc_params(); end
if isempty(w0), w0 = p.w0; end
if isempty(q0), q0 = p.q0; end
q0 = q0(:)/norm(q0);
x0 = [w0(:); q0];

switch lower(o.controller)
    case 'pakmpc'
        c = koopman_mpc(quat_koopman(p.nphi, p.jc0, p.ts), p, struct('adapt', true));
    case 'nonadaptive'
        c = koopman_mpc(quat_koopman(p.nphi, p.jc0, p.ts), p);
    case 'kmae'
        c = kmae_akmpc(p);
    case 'inmpc'
        c = incremental_nmpc(p);
end

plant = struct('tend', o.tend, 'jc', o.jc, 'target', o.target, 'noise', o.noise, ...
    'seed', o.seed, 'target_ref', 'identity');      % [guess] as in fig6_11_table3_adaptive
if o.dist, plant.dist = p.dist; end
out = simulate_attitude(c, x0, p, plant);

s = out.t >= out.t(end) - 10;
out.settle = settling_time(out.t, out.x(5:7, :));
out.ss_qv = mean(vecnorm(out.x(5:7, s)));
out.ss_w = mean(vecnorm(out.x(1:3, s)));
out.effort = sum(vecnorm(out.tau).^2)*p.ts;
fprintf('%s: settled (5%%) at %.1f s, steady-state ||q_v|| = %.2e, ||w|| = %.2e, effort = %.3g, step %.1f ms\n', ...
    o.controller, out.settle, out.ss_qv, out.ss_w, out.effort, 1e3*mean(out.steptime));

if o.plot
    adapt = strcmpi(o.controller, 'pakmpc');
    figure('Name', ['pakmpc: ' o.controller]);
    tl = tiledlayout(3 + 2*adapt, 1, 'TileSpacing', 'compact');
    nexttile; plot(out.t, out.x(5:7, :)); grid on; ylabel('q_v'); legend('q_1', 'q_2', 'q_3');
    nexttile; plot(out.t, out.x(1:3, :)); grid on; ylabel('\omega (rad/s)');
    legend('\omega_1', '\omega_2', '\omega_3');
    nexttile; stairs(out.t(1:end-1), out.tau'); grid on; ylabel('\tau (N\cdotm)');
    if adapt
        nexttile; plot(out.t(1:end-1), out.mu); grid on; ylabel('\mu_k');
        nexttile; plot(out.t(1:end-1), out.db); grid on; ylabel('||\DeltaB||_2');
    end
    xlabel(tl, 'Time (s)');
end
end
