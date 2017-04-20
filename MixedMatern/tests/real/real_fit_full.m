% run on server
parpool(16)
addpath(genpath('/home/minjay/div_curl'))

% full bivariate Matern model with a1=a2=a12

load('wind.mat')

% initial computation
[h_mat, r, P_cell, Q_cell, A_cell] = init_comp(x, y, z, n, theta, phi);

% negative log-likelihood function
negloglik1 = @(beta_all) negloglik_full(beta_all, h_mat, r, P_cell, Q_cell, A_cell, samples);

beta_init = [0.029113 0.054503 0.281047 1.757812 2.034312 9.47189 0.210496 0.196141];
beta_init = [beta_init(1:5) (beta_init(4)+beta_init(5))/2 beta_init(6:end)];

lb = [0 0 -1 1 1 1 0 0 0];
ub = [Inf Inf 1 5 5 5 Inf Inf Inf];

% fit the model
[beta_hat, f_min] = Matern_fit(negloglik1, beta_init, lb, ub, @mycon_full, true);

delete(gcp)
