% simulate the data for examining the prediction performance

clear

% HEALPix grid
[theta, phi, n] = HEALPix_sampling(3);

% convert coordinates
[x, y, z] = trans_coord(theta, phi);

% initial computation
[h_mat, r, P_cell, Q_cell, A_cell] = init_comp(x, y, z, n, theta, phi);

p = 2;
% set the same seed
rng('default')

% specify parameters
% [sigma1, sigma2, rho12, nu1, nu2, a, tau1, tau2]
beta_all = [1 1 0.5 3 4 2 0.1 0.1];

samples = mvnrnd(zeros(p*n, 1), eye(p*n))';
[u, v] = sim_mix(n, beta_all, samples, h_mat, r, P_cell, Q_cell, A_cell);

filename = 'sim_data_mix.mat';
save(filename, 'u', 'v', 'theta', 'phi', 'n', 'x', 'y', 'z')
