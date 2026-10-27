function c = kmae_akmpc(p, opts)
%KMAE_AKMPC  KMAE-AKMPC baseline of Sec. V-C, after Chen et al. [31]: linear MPC on an RBF-EDMD
%   Koopman model whose A and B are updated online by recursive least squares (RLS).
%   c = kmae_akmpc(p, opts); then each sample [tau, c, info] = c.step(c, x), as koopman_mpc.
%   opts fields (all optional):
%     lam_f     RLS forgetting factor, default 0.995                                     [guess]
%     interp    interpolated initial state of koopman_mpc, default false                 [guess]
%
%   [paper] "a linear MPC based on a nominal model, which is initially established using the EDMD
%   method and then continuously updated online via a Recursive Least Squares (RLS) algorithm".
%   [guess] The nominal model is rbf_edmd (50 thin-plate RBFs) trained on the nominal inertia
%   Jc0 = Js, so it starts with the same model error as PAKMPC. RLS continues the offline least
%   squares fit: Theta = [A B] starts at the EDMD solution and P0 = inv(Gram) of the training data,
%   and each sample the pair (z_(k-1), tau_(k-1)) -> z_k updates Theta. The MPC (cost, weights,
%   horizon, terminal weight, input bounds) is koopman_mpc with the updated (A, B).

if nargin < 2, opts = struct(); end
d = struct('lam_f', 0.995, 'interp', false);
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(opts, f{i}) || isempty(opts.(f{i})), opts.(f{i}) = d.(f{i}); end
end
m = rbf_edmd(p, 50, 1, p.jc0);
m.bconst0 = m.bconst;                       % offline B, for info.db
c.inner = koopman_mpc(m, p, struct('interp', opts.interp));
c.theta = [m.a m.bconst];
c.pr = inv(m.gram + 1e-9*eye(size(m.gram)));
c.lam_f = opts.lam_f;
c.prev = [];
c.step = @step;
end

function [tau, c, info] = step(c, x)
z = c.inner.m.lift(x);
if ~isempty(c.prev)
    % RLS on z_k = Theta [z_(k-1); tau_(k-1)]
    r = c.prev;
    k = c.pr*r/(c.lam_f + r'*c.pr*r);
    c.theta = c.theta + (z - c.theta*r)*k';
    c.pr = (c.pr - k*(r'*c.pr))/c.lam_f;
    c.pr = (c.pr + c.pr')/2;
    c.inner = set_model(c.inner, c.theta);
end
[tau, c.inner, info] = c.inner.step(c.inner, x);
c.prev = [z; tau];
info.db = norm(c.theta(:, end-2:end) - c.inner.m.bconst0);
end

function ci = set_model(ci, theta)
% swap (A, B) in the koopman_mpc struct and rebuild its prediction matrices
m = ci.m; n = m.n; N = ci.N;
m.a = theta(:, 1:n);
m.bconst = theta(:, n+1:end);
bc = m.bconst;
m.b = @(w, q) bc;
ak = eye(n);
for i = 0:N
    ci.ab(i*n + (1:n), :) = ak;
    ak = m.a*ak;
end
ci.m = m;
end
