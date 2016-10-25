% MLE of the Mixed Matern model
% number of grid points=768

clear

savefile = 'MLE_healpix.mat';

% run on server
% 24 is the maximum on hannan
parpool(24)
addpath(genpath('/home/minjay/div_curl'))
addpath(genpath('/home/minjay/MEALPix'))

% HEALPix grid
[theta, phi, n] = HEALPix_sampling(3);

% convert coordinates
[x, y, z] = trans_coord(theta, phi);

% initial computation
[h_mat, r, P_cell, Q_cell, A_cell] = init_comp(x, y, z, n, theta, phi);

p = 2;
% set the same seed
rng('default')

N = 100;
lb = [0 0 -1 1 1 0 0 0];
ub = [Inf Inf 1 10 10 Inf Inf Inf];

% specify parameters
% [sigma1, sigma2, rho12, nu1, nu2, a, tau1, tau2]
beta_all = [1 1 0.5 3 4 2 0.1 0.1];
rec_beta_hat = zeros(N, length(beta_all));

parfor rep = 1:N
    samples = mvnrnd(zeros(p*n, 1), eye(p*n))';
    [u, v] = sim_mix(n, beta_all, samples, h_mat, r, P_cell, Q_cell, A_cell);
    samples = reshape([u v]', 1, p*n);

    % negative log-likelihood function
    negloglik1 = @(beta_all) negloglik(beta_all, h_mat, r, P_cell, Q_cell, A_cell, samples);

    % rand search
    beta_init = rand_search(negloglik1, 100, [0 0 -1 1 1 0 0 0], [10 10 1 10 10 10 10 10], true, false);

    % fit the model
    [beta_hat, f_min] = Matern_fit(negloglik1, beta_init, lb, ub, @mycon, false);
    
    rec_beta_hat(rep, :) = beta_hat;
end

save(savefile, 'rec_beta_hat');

% run on server
delete(gcp)
