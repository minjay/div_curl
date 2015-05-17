function [] = sim4_bootstrap(num)

% run on server
parpool(8)
addpath(genpath('/home/minjay/div_curl'))

load('MLE2.mat')

savefile = ['boot_sim', num2str(num), '.mat'];

p = 2;
% number of bootstrapping
B = 200;

n_lat = 15;
n_lon = 30;

% regular grid
[theta, phi, n] = regular_sampling(n_lat, n_lon);

% convert coordinates
[x, y, z] = trans_coord(theta, phi);

% initial computation
[h_mat, r, P_cell, Q_cell, A_cell] = init_comp(x, y, z, n, theta, phi);

beta_all = rec_beta_hat(num, :);
beta = beta_all(1:6);
tau1 = beta_all(7);
tau2 = beta_all(8);

[coef, bessel] = get_coef_bessel(beta, r);

% get cov mat
cov_mat = get_cov(h_mat, r, P_cell, Q_cell, A_cell, @Matern_mix,...
    beta, coef, bessel)+diag(kron(ones(1, n), [tau1^2, tau2^2]));

samples_all = mvnrnd(zeros(p*n, 1), cov_mat, B);

beta_init = beta_all;
lb = [0 0 -1 1 1 0 0 0];
ub = [Inf Inf 1 5 5 Inf Inf Inf];

rec_beta_hat = zeros(B, 8);

% fit the model
parfor rep = 1:B
    
    samples = samples_all(rep, :);
    
    % negative log-likelihood function
    negloglik1 = @(beta_all) negloglik_fast(beta_all, h_mat, r, P_cell, Q_cell, A_cell, samples, n_lat, n_lon);
    
    [beta_hat, f_min] = Matern_fit(negloglik1, beta_init, lb, ub, @mycon, false);
    rec_beta_hat(rep, :) = beta_hat;
    
end

save(savefile, 'rec_beta_hat');

% run on server
delete(gcp)

end
