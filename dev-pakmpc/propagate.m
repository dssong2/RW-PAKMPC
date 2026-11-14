p = params();
t = 0:p.dt:100;
x = zeros(7, numel(t));
x(:, 1) = [p.w0; p.q0];
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

disturb = false;
disturbance = @(t) [0.1*sin(0.075*t); -0.2*sin(0.1*t); 0.15*sin(0.1*t)];

for k = 1:numel(t) - 1
    % Define input, replaced with controller later
    tau = [0; 0; 0];

    % odefun selections: no disturbance, with disturbance
    odefun = @(s, y) attitude_ode(y, tau, p.Jc);
    % Disturbance variation
    odefun_dist = @(s, y) attitude_ode(y, tau + disturb*disturbance(t(k) + s), p.Jc);

    % Solve and propagate with ode45
    [~, xx] = ode45(odefun_dist, [0 p.dt], x(:, k), opts);
    x(:, k+1) = xx(end, :)';

    % Renormalize quaternions
    x(4:7, k + 1) = x(4:7, k + 1)/norm(x(4:7, k + 1));
end

% Claude
R = @(q) (q(1)^2 - q(2:4)'*q(2:4))*eye(3) + 2*(q(2:4)*q(2:4)') ...
    + 2*q(1)*[0 -q(4) q(3); q(4) 0 -q(2); -q(3) q(2) 0];
H = zeros(3, numel(t)); E = zeros(1, numel(t));
for k = 1:numel(t)
    w = x(1:3, k);
    H(:, k) = R(x(4:7, k))*p.Jc*w;
    E(k) = 0.5*w'*p.Jc*w;
end
fprintf('max |H - H0|/|H0| = %.2e\n', max(vecnorm(H - H(:, 1)))/norm(H(:, 1)));
fprintf('max |E - E0|/E0   = %.2e\n', max(abs(E - E(1)))/E(1));