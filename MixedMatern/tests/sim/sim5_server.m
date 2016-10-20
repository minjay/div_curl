% MLE of the Mixed Matern model
% In this script, we consider the case where covariates are included in the
% model.
% The covariates we included are sin(longitude) and cos(latitude).
% grid=20*40

clear

savefile = 'MLE_covariate.mat';

% run on server
parpool(10)
addpath(genpath('/home/minjay/div_curl'))
addpath(genpath('/home/minjay/MEALPix'))

n_lat = 20;
n_lon = 40;

% regular grid
[theta, phi, n] = regular_sampling(n_lat, n_lon);

% convert coordinates
[x, y, z] = trans_coord(theta, phi);

% initial computation
[h_mat, r, P_cell, Q_cell, A_cell] = init_comp(x, y, z, n, theta, phi);

p = 2;
% set the same seed
rng('default')

N = 500;
lb = [0 0 -1 1 1 0 0 0 -Inf, -Inf, -Inf, -Inf, -Inf, -Inf];
ub = [Inf Inf 1 5 5 Inf Inf Inf, Inf, Inf, Inf, Inf, Inf, Inf];

% specify parameters
% [sigma1, sigma2, rho12, nu1, nu2, a, tau1, tau2, c10, c11, c12, c20, c21, c22]
beta_all = [1 1 0.5 3 4 2 0.1 0.1, 5, 5, -1, 10, -15, 3];
beta_partial = beta_all(1:8);
rec_beta_hat = zeros(N, length(beta_all));
samples_all = mvnrnd(zeros(p*n, 1), eye(p*n), N);
% specify mean component
m_u = beta_all(9) + beta_all(10) * theta + beta_all(11) * theta.^2;
m_v = beta_all(12) + beta_all(13) * theta + beta_all(14) * theta.^2;

parfor rep = 1:N
    samples = samples_all(rep, :)';
    [u, v] = sim_mix(n, beta_partial, samples, h_mat, r, P_cell, Q_cell, A_cell);
    
    % add mean component
    u = u + m_u;
    v = v + m_v;
    
    % use linear regression to get initial values of coefficient c's
    X = [ones(n, 1) theta theta.^2];
    coef_u = (X' * X) \ (X' * u);
    coef_v = (X' * X) \ (X' * v);
    
    samples = reshape([u v]', 1, p*n);

    % negative log-likelihood function
    negloglik1 = @(beta_all) negloglik_fast_covariate(beta_all, h_mat, r, P_cell, Q_cell, A_cell, samples, theta, phi, n_lat, n_lon);

    % rand search
    beta_init = rand_search(negloglik1, 100, [0 0 -1 1 1 5 1e-3 1e-3 coef_u' coef_v'], [5 5 1 5 5 5 1e-3 1e-3 coef_u' coef_v'], true, false);

    % fit the model
    [beta_hat, f_min] = Matern_fit(negloglik1, beta_init, lb, ub, @mycon, false);
    
    rec_beta_hat(rep, :) = beta_hat;
end

save(savefile, 'rec_beta_hat');

% run on server
delete(gcp)
