% MLE of the Mixed Matern model
% In this script, we consider the case where covariates are included in the
% model.
% The covariates we included are sin(longitude) and cos(latitude).
% grid=10*20

clear

savefile = 'MLE_covariate_good_init.mat';

% run on server
parpool(16)
addpath(genpath('/home/minjay/div_curl'))
addpath(genpath('/home/minjay/MEALPix'))

n_lat = 10;
n_lon = 20;

% regular grid
[theta, phi, n] = regular_sampling(n_lat, n_lon);

% convert coordinates
[x, y, z] = trans_coord(theta, phi);

% initial computation
[h_mat, r, P_cell, Q_cell, A_cell] = init_comp(x, y, z, n, theta, phi);

p = 2;
% set the same seed
rng('default')

N = 100;
lb = [0 0 -1 1 1 0 0 0 -Inf, -Inf, -Inf, -Inf, -Inf, -Inf];
ub = [Inf Inf 1 5 5 Inf Inf Inf, Inf, Inf, Inf, Inf, Inf, Inf];

% specify parameters
% [sigma1, sigma2, rho12, nu1, nu2, a, tau1, tau2, c10, c11, c12, c20, c21, c22]
beta_all = [1 1 0.5 3 4 2 0.1 0.1, 5, 5, -1, 10, -15, 3];
beta_init = beta_all;
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
    
    samples = reshape([u v]', 1, p*n);

    % negative log-likelihood function
    negloglik1 = @(beta_all) negloglik_fast_covariate(beta_all, h_mat, r, P_cell, Q_cell, A_cell, samples, theta, phi, n_lat, n_lon);

    % fit the model
    [beta_hat, f_min] = Matern_fit_long(negloglik1, beta_init, lb, ub, @mycon, false);
    
    rec_beta_hat(rep, :) = beta_hat;
end

save(savefile, 'rec_beta_hat');

% run on server
delete(gcp)
