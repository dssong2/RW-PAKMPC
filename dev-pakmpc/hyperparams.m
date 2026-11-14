%% Table 1 Practical Hyper-Parameter Selection
function hp = hyperparams()
hp.N_phi = 5; % Koopman Truncation Order
hp.Hw = 10; % data window length
hp.N = 15; % prediction horizon
hp.lambda = 0.965; % forgetting factor
hp.rho = 1e-7; % learning rate
hp.gamma = 1e-3; % regularization weight
hp.beta1 = 0.9; % Adam first moment
hp.beta2 = 0.999; % Adam second moment
hp.epsilon = 1e-8; % Adam small constant
end