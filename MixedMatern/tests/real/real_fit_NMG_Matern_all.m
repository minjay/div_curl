% run on server
parpool(8)
addpath(genpath('/home/minjay/div_curl'))

load('wind.mat')

% initial computation
[r, h0_cell] = init_comp_NMG(n, theta, phi, x, y, z);

param_BM = [0.21810846 0.203180249 0.058196367 1.23929461 1.1322247 0.184217146 0.156812964 -0.07977680];
beta_all(1:2) = param_BM(6:7);
beta_all(3) = param_BM(8);
beta_all(4:5) = param_BM(4:5);
beta_all(6) = param_BM(3);
beta_all(7:8) = param_BM(1:2);
sigma1 = sqrt(beta_all(1));
sigma2 = sqrt(beta_all(2));
rho12 = beta_all(3);
nu1 = beta_all(4);
nu2 = beta_all(5);
a = 1/beta_all(6);
tau1 = beta_all(7);
tau2 = beta_all(8);

% negative log-likelihood function
negloglik1 = @(beta_all) negloglik_NMG_Matern_all(beta_all, r, samples, h0_cell);

% [a1 a2 b1 b2 nu a sigma1 sigma2 w1 w2 tau1 tau2]
beta_init = [0.006925 -0.003031 0.000483 0.005205 3.305484 a sigma1 sigma2 nu1 nu2 tau1 tau2];
negloglik1(beta_init)

% to avoid identifiability problem, set a1>0
lb = [0   -Inf -Inf -Inf 1 0   0   0   1 1 0   0];
ub = [Inf Inf  Inf  Inf  5 Inf Inf Inf 5 5 Inf Inf];

% fit the model
[beta_hat, f_min] = Matern_fit(negloglik1, beta_init, lb, ub, [], true);

delete(gcp)
