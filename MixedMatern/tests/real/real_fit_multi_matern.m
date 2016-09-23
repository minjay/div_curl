% run on server
parpool(8)
addpath(genpath('/home/minjay/div_curl'))

load('wind.mat')
load('samples_multi_matern.mat')

% initial computation
[h_mat, r, P_cell, Q_cell, A_cell] = init_comp(x, y, z, n, theta, phi);

% negative log-likelihood function
negloglik1 = @(beta_all) negloglik_multi_matern(beta_all, r, samples_all(1:108, :));

beta_init = [0.21810846 0.203180249 0.058196367 1.23929461 1.1322247 0.184217146 0.156812964 -0.07977680];

beta_all(1:2) = beta_init(6:7);
beta_all(3) = beta_init(8);
beta_all(4:5) = beta_init(4:5);
beta_all(6) = beta_init(3);
beta_all(7:8) = beta_init(1:2);
sigma1 = sqrt(beta_all(1));
sigma2 = sqrt(beta_all(2));
rho12 = beta_all(3);
nu1 = beta_all(4);
nu2 = beta_all(5);
a = 1/beta_all(6);
tau1 = beta_all(7);
tau2 = beta_all(8);

beta_init = [sigma1 sigma2 rho12 nu1 nu2 a tau1 tau2];

lb = [0 0 -1 0 0 0 0 0];
ub = [Inf Inf 1 5 5 Inf Inf Inf];

% fit the model
[beta_hat, f_min] = Matern_fit(negloglik1, beta_init, lb, ub, @mycon, true);

delete(gcp)
