% consider the data with 90% EOF subtracted

% run on server
parpool(8)
addpath(genpath('/home/minjay/div_curl'))

load('wind.mat')

% initial computation
[h_mat, r, P_cell, Q_cell, A_cell] = init_comp(x, y, z, n, theta, phi);

% negative log-likelihood function
negloglik1 = @(beta_all) negloglik(beta_all, h_mat, r, P_cell, Q_cell, A_cell, samples);

% rand search
beta_init = rand_search(negloglik1, 100, [0 0 -1 1 1 10 0 0], [10 10 1 5 5 10 10 10], true, true);

lb = [0 0 -1 1 1 0 0 0];
ub = [Inf Inf 1 5 5 Inf Inf Inf];

% fit the model
[beta_hat, f_min] = Matern_fit(negloglik1, beta_init, lb, ub, @mycon, true);

delete(gcp)
