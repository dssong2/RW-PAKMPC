function sim_out = run_sim_scenario(target_scale, dist_scale, mode, mdl_name)
    % RUN_SIM_SCENARIO Master execution wrapper for Simscape evaluations.
    %
    % Inputs:
    %   target_scale - Scalar multiplier for target mass & inertia (e.g., 0.5 to 1.5)
    %   dist_scale   - Scalar multiplier for external disturbance & target torque
    %   mode         - 'rw' (Tetrahedral Reaction Wheels) or 'ideal' (Direct Torque)
    %   mdl_name     - Name of Simulink model (default: 'pakmpc_rw_combined')

    if nargin < 4 || isempty(mdl_name)
        mdl_name = 'pakmpc_rw_combined';
    end
    if nargin < 3 || isempty(mode)
        mode = 'rw';
    end

    % 1. Ensure baseline configuration exists in Base Workspace
    if evalin('base', '~exist(''ms'', ''var'')')
        evalin('base', 'config;');
    end

    % Pull nominal parameters from base workspace
    mt_nom      = evalin('base', 'mt');
    jt_nom      = evalin('base', 'jt');
    jw_val      = evalin('base', 'jw');
    wt_val      = evalin('base', 'wt');
    ts_val      = evalin('base', 'ts');
    tend_val    = evalin('base', 'tend');
    noise_on    = evalin('base', 'noise_on');
    target_on   = evalin('base', 'target_on');
    dist_on     = evalin('base', 'dist_on');

    % 2. Calculate scaled target physical properties
    mt_actual = mt_nom * target_scale;
    jt_actual = jt_nom * target_scale;

    % Load model silently if not already loaded
    if ~bdIsLoaded(mdl_name)
        load_system(mdl_name);
    end

    % 3. Update Target Inertia block in Simscape
    set_param([mdl_name '/Target'], 'Mass', num2str(mt_actual, 17));
    set_param([mdl_name '/Target'], 'MomentsOfInertia', mat2str(diag(jt_actual)', 17));
    set_param([mdl_name '/Target'], 'ProductsOfInertia', ...
        mat2str([jt_actual(2,3), jt_actual(3,1), jt_actual(1,2)], 17));

    % 4. Update Disturbance & Target Torque Gains for sweep
    set_param([mdl_name '/Target On'], 'Gain', num2str(target_on * dist_scale, 17));
    set_param([mdl_name '/Disturbance On'], 'Gain', num2str(dist_on * dist_scale, 17));

    % 5. Configure Actuation Mode ('rw' vs 'ideal')
    if strcmpi(mode, 'ideal')
        set_param([mdl_name '/Ideal Direct Mode'], 'Gain', '1');
        set_param([mdl_name '/RW Mode'], 'Gain', '0');
        set_param([mdl_name '/RW Dist On'], 'Gain', '0');
    else
        set_param([mdl_name '/Ideal Direct Mode'], 'Gain', '0');
        set_param([mdl_name '/RW Mode'], 'Gain', '1');
        set_param([mdl_name '/RW Dist On'], 'Gain', num2str(noise_on));
    end

    % 6. Run simulation and time execution
    t_wall_start = tic;
    out = sim(mdl_name, 'StopTime', num2str(tend_val), 'ReturnWorkspaceOutputs', 'on');
    wall_time = toc(t_wall_start);

    % 7. Extract and format logged signals
    sim_out.mode         = mode;
    sim_out.target_scale = target_scale;
    sim_out.dist_scale   = dist_scale;
    sim_out.t            = out.q_log.time(:);
    num_t                = length(sim_out.t);

    sim_out.q = format_sig(out.q_log, num_t, 4);
    sim_out.w = format_sig(out.w_log, num_t, 3);
    sim_out.p = format_sig(out.p_log, num_t, 4);

    % Discrete signals (sampled at ts)
    sim_out.tm     = out.mm_log.time(:);
    num_tm         = length(sim_out.tm);
    sim_out.M_m_d  = format_sig(out.mm_log, num_tm, 4);
    sim_out.u_d    = format_sig(out.u_log, num_tm, 3);

    % Resample discrete commands onto continuous time t (Zero-Order Hold)
    sim_out.M_m = interp1(sim_out.tm, sim_out.M_m_d, sim_out.t, 'previous', 'extrap');
    sim_out.u   = interp1(sim_out.tm, sim_out.u_d,   sim_out.t, 'previous', 'extrap');

    % Derived Reaction Wheel telemetry
    sim_out.hw         = jw_val * sim_out.p;                  % Wheel angular momentum (N*m*s)
    sim_out.tau_rw_net = (-wt_val * sim_out.M_m.')';          % Net 3-axis RW body torque (N*m)

    % Solver timing metrics
    num_discrete_steps     = floor(tend_val / ts_val);
    sim_out.wall_time      = wall_time;
    sim_out.num_steps      = num_discrete_steps;
    sim_out.avg_step_time  = wall_time / num_discrete_steps;  % Average step runtime (s)
end

function x = format_sig(log_struct, expected_rows, expected_cols)
    x = squeeze(log_struct.signals.values);
    if size(x, 1) ~= expected_rows && size(x, 2) == expected_rows
        x = x.';
    end
    if size(x, 2) ~= expected_cols && size(x, 1) == expected_cols
        x = x.';
    end
end