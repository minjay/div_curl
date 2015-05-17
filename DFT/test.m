% a script to test:
% (1) whether negloglik_fast is correct
% (2) how fast negloglik_fast runs

clear

n_lat = 25;
n_lon = 50;

% regular grid
[theta, phi, n] = regular_sampling(n_lat, n_lon);

% convert coordinates
[x, y, z] = trans_coord(theta, phi);

% initial computation
[h_mat, r, P_cell, Q_cell, A_cell] = init_comp(x, y, z, n, theta, phi);

% generate samples
p = 2;
% set the same seed
rng('default')
samples = mvnrnd(zeros(p*n, 1), eye(p*n))';

% specify parameters
% [sigma1, sigma2, rho12, nu1, nu2, a, tau1, tau2]
beta_all = [1 1 0.5 3 4 2 0.1 0.1];
[u, v] = sim_mix(n, beta_all, samples, h_mat, r, P_cell, Q_cell, A_cell);
samples = reshape([u v]', 1, p*n);

negloglik1 = @(beta_all) negloglik(beta_all, h_mat, r, P_cell, Q_cell, A_cell, samples);
negloglik2 = @(beta_all) negloglik_fast(beta_all, h_mat, r, P_cell, Q_cell, A_cell, samples, n_lat, n_lon);

tic
negloglik1(beta_all)
toc
tic
negloglik2(beta_all)
toc
