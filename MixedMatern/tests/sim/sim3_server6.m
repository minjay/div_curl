% MLE of the Mixed Matern model
% record the computational time of the case with 5000 data points
% no parallel computation
% the only "true" is for LHS sampling

clear

% run on server
addpath(genpath('/home/minjay/div_curl'))
addpath(genpath('/home/minjay/MEALPix'))

n_lat = 50;
n_lon = 100;

% regular grid
[theta, phi, n] = regular_sampling(n_lat, n_lon);

% convert coordinates
[x, y, z] = trans_coord(theta, phi);

% initial computation
[h_mat, r, P_cell, Q_cell, A_cell] = init_comp(x, y, z, n, theta, phi);

p = 2;
% set the same seed
rng('default')

N = 10;
lb = [0 0 -1 1 1 0 0 0];
ub = [Inf Inf 1 5 5 Inf Inf Inf];

% specify parameters
% [sigma1, sigma2, rho12, nu1, nu2, a, tau1, tau2]
beta_all = [1 1 0.5 3 4 2 0.1 0.1];
samples_all = mvnrnd(zeros(p*n, 1), eye(p*n), N);
time_all = zeros(N, 1);

for rep = 1:N
    samples = samples_all(rep, :)';
    [u, v] = sim_mix(n, beta_all, samples, h_mat, r, P_cell, Q_cell, A_cell);
    samples = reshape([u v]', 1, p*n);

    % negative log-likelihood function
    negloglik1 = @(beta_all) negloglik_fast(beta_all, h_mat, r, P_cell, Q_cell, A_cell, samples, n_lat, n_lon);

    % rand search
    beta_init = rand_search(negloglik1, 100, [0 0 -1 1 1 5 1e-3 1e-3], [5 5 1 5 5 5 1e-3 1e-3], true, false);
    
    tic
    % fit the model
    [beta_hat, f_min] = Matern_fit(negloglik1, beta_init, lb, ub, @mycon, false); 
    time_all(rep) = toc;
end

save('time_5000.mat', 'time_all');
