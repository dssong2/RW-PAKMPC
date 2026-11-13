function rw = rw_params(preset)
%RW_PARAMS  All reaction wheel parameters for the tetrahedral 4-wheel array.
%   rw = rw_params                   default wheel, 'realistic_hr16'
%   rw = rw_params('realistic_hr16') Honeywell HR16-50: 0.2 N m, 50 N m s, 9.0 kg [datasheet]
%   rw = rw_params('hr16_100')       Honeywell HR16-100: 0.2 N m, 100 N m s, 12.0 kg [datasheet]
%   rw = rw_params('high_torque')    15 N m, 100 N m s: sized to the paper's 20 N m authority, not a real wheel
%
%   Datasheet: Honeywell, "Constellation Series Reaction Wheels", N61-0045-001-001, Dec 2003.
%   HR16: 50/75/100 N m s at 6000 rpm, reaction torque 0.1 to 0.2 N m nominal (up to 0.4 extended
%   at 3000 rpm), unit mass 9/10.4/12 kg.
%
%   Used by build_rw_array (block defaults and checks) and build_combined_multibody (rw_preset).
if nargin < 1, preset = 'realistic_hr16'; end

% ---------------- geometry ----------------
beta = acos(1/sqrt(3));
rw.wt = [sin(beta) -sin(beta) 0          0;
         cos(beta)  cos(beta) cos(beta)  cos(beta);
         0          0         sin(beta) -sin(beta)];   % spin axes in the body frame
rw.wp = pinv(rw.wt);                    % = 0.75*W_T'
rw.n = [1; 1; -1; -1]/2;                % unit null-space vector, W_T*n = 0
rw.s = 0.32;                            % wheel spacing, m [peer's cluster radius]
rw.d = rw.s*(rw.wt - mean(rw.wt, 2));   % rotor positions in B, sum to zero so the CoM stays at B
rw.rg = 0.16; rw.lg = 0.06;             % graphics only: rotor radius and thickness, m

% ---------------- wheel sizing ----------------
rpm2rad = 2*pi/60;
switch preset
    case 'realistic_hr16'               % HR16-50 [datasheet]
        rw.hmax   = 50;                 % N m s at 6000 rpm
        rw.jw     = rw.hmax/(6000*rpm2rad);   % rotor axial inertia = h_max / max speed, 0.0796 kg m^2
        rw.mw     = 9.0;                % whole-unit mass, kg, lumped at the rotor [lumping assumed]
        rw.taumax = 0.2;                % top of the nominal range, held constant over speed [constant assumed]
    case 'hr16_100'                     % HR16-100 [datasheet]
        rw.hmax   = 100;
        rw.jw     = rw.hmax/(6000*rpm2rad);
        rw.mw     = 12.0;
        rw.taumax = 0.2;
    case 'high_torque'                  % peer's upper bound, not a real wheel
        rw.jw     = 0.15;
        rw.mw     = 8.0;
        rw.taumax = 15.0;
        rw.hmax   = 100;
    otherwise
        error('rw_params: unknown preset ''%s''.', preset);
end
rw.jt = rw.jw/2;                        % rotor transverse inertia, thin disk [assumed]

% ---------------- friction and torque error ----------------
rw.cth = 1e-4;                          % viscous bearing drag, N m s/rad [assumed]
rw.tc = 0;                              % Coulomb friction, N m [assumed, off]
rw.pc = 0.1;                            % tanh smoothing speed for Coulomb friction, rad/s
rw.sigma = 1e-4;                        % torque noise std, N m [assumed]
rw.b = 2e-4*[1 -1 0.5 -0.5];            % constant torque bias, N m [assumed]
rw.sigd = 1e-5;                         % drift random walk, N m/sqrt(s) [assumed]
rw.dist = 1;                            % 0 turns noise, bias and drift off
rw.seed = 0;                            % noise seeds seed+1:seed+4, drift seeds seed+11:seed+14
rw.ts = 0.1;                            % controller sample time, s [paper]

% ---------------- mass properties added to the spacecraft ----------------
rw.mass = 4*rw.mw;
rw.ilock = zeros(3);                    % locked-rotor inertia about B
for i = 1:4
    w = rw.wt(:,i); d = rw.d(:,i);
    rw.ilock = rw.ilock + rw.jt*(eye(3) - w*w') + rw.jw*(w*w') + rw.mw*((d'*d)*eye(3) - d*d');
end
end
