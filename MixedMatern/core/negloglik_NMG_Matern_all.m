% run on server
parpool(16)
addpath(genpath('/home/minjay/div_curl'))

load('wind.mat')

% initial computation
[r, h0_cell] = init_comp_NMG(n, theta, phi, x, y, z);

% negative log-likelihood function
negloglik1 = @(beta_all) negloglik_NMG_Matern(beta_all, r, samples, h0_cell);

% [a1 a2 b1 b2 nu a sigma1 sigma2 w1 w2 tau1 tau2]
beta_init = [0.090899 -0.056781 -0.00133 0.047889 2.735421 17.183203 0.429205 0.395996 1.239295 1.132225 0.218108 0.20318];
negloglik1(beta_init)

% to avoid identifiability problem, set a1>0
lb = [0   -Inf -Inf -Inf 1 0   0   0   1 1 0   0];
ub = [Inf Inf  Inf  Inf  5 Inf Inf Inf 5 5 Inf Inf];

% fit the model
[beta_hat, f_min] = Matern_fit(negloglik1, beta_init, lb, ub, [], true);

delete(gcp)
