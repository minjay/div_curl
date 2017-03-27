% run on server
parpool(8)
addpath(genpath('/home/minjay/div_curl'))

load('sim_data_mix.mat')
load('param_kriging_sim.mat')

% initial computation
[r, h0_cell] = init_comp_NMG(n, theta, phi, x, y, z);

rep = 1;
beta_all(1:2) = param_BM(rep, 6:7);
beta_all(3) = param_BM(rep, 8);
beta_all(4:5) = param_BM(rep, 4:5);
beta_all(6) = param_BM(rep, 3);
beta_all(7:8) = param_BM(rep, 1:2);
sigma1 = sqrt(beta_all(1));
sigma2 = sqrt(beta_all(2));
rho12 = beta_all(3);
nu1 = beta_all(4);
nu2 = beta_all(5);
a = 1/beta_all(6);
tau1 = beta_all(7);
tau2 = beta_all(8);
cov_mat_Matern = get_cov_Matern_pars(r, sigma1, sigma2, 0, nu1, nu2, a);

% negative log-likelihood function
negloglik1 = @(beta_all) negloglik_NMG_Matern(beta_all, r, samples, h0_cell, cov_mat_Matern, a, tau1, tau2);

% [a1 a2 b1 b2 nu a sigma1 sigma2 w1 w2 tau1 tau2]
beta_init = [0 0 0 0 2+1];
negloglik1(beta_init)

% to avoid identifiability problem, set a1>0
lb = [0   -Inf -Inf -Inf 1];
ub = [Inf Inf  Inf  Inf  5];

% fit the model
[beta_hat, f_min] = Matern_fit(negloglik1, beta_init, lb, ub, [], true);

delete(gcp)