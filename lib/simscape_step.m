function tau = simscape_step(wm, qm)
%SIMSCAPE_STEP  Bridge from the Simscape model's Controller block to a koopman_mpc-style controller.
%   simscape_step('init', c)   stores controller c (anything with [tau, c, info] = c.step(c, x))
%   tau = simscape_step(wm, qm)   one controller sample on the measured w and q, called through
%                                 coder.extrinsic from the Controller MATLAB Function block
%   lg = simscape_step('log')  returns the per-sample log: applied tau, mu, ||Delta B||, wall time

persistent c lg
if ischar(wm)
    switch wm
        case 'init'
            c = qm;
            lg = struct('tau', zeros(3, 0), 'mu', [], 'db', [], 'steptime', []);
        case 'log'
            tau = lg;
    end
    return
end
t0 = tic;
[tau, c, info] = c.step(c, [wm(:); qm(:)]);
lg.steptime(end + 1) = toc(t0);
lg.tau(:, end + 1) = tau;
lg.mu(end + 1) = info.mu;
lg.db(end + 1) = info.db;
end
