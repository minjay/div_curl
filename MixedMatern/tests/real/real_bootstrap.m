% compute the standard errors by bootstrapping
clear

% run on server
parpool(16)
addpath(genpath('/home/minjay/div_curl'))

load('wind.mat')

savefile = 'boot1.mat';

T = 108;
p = 2;
% number of bootstrapping
B = 200;

% initial computation
[h_mat, r, P_cell, Q_cell, A_cell] = init_comp(x, y, z, n, theta, phi);

beta_all = [0.029113 0.054503 0.281047 1.757812 2.034312 9.47189 0.210496 0.196141];
beta = beta_all(1:6);
tau1 = beta_all(7);
tau2 = beta_all(8);

[coef, bessel] = get_coef_bessel(beta, r);

% get cov mat
cov_mat = get_cov(h_mat, r, P_cell, Q_cell, A_cell, @Matern_mix,...
    beta, coef, bessel)+diag(kron(ones(1, n), [tau1^2, tau2^2]));

samples_all = mvnrnd(zeros(p*n, 1), cov_mat, B*T);

samples_all_cell = cell(B, 1);
for rep = 1:B
    samples_all_cell{rep} = samples_all((rep-1)*T+1:rep*T, :);
end

beta_init = beta_all;
lb = [0 0 -1 1 1 0 0 0];
ub = [Inf Inf 1 5 5 Inf Inf Inf];

rec_beta_hat = zeros(B, 8);

% fit the model
parfor rep = 1:B
    
    samples = samples_all_cell{rep};
    
    % negative log-likelihood function
    negloglik1 = @(beta_all) negloglik(beta_all, h_mat, r, P_cell, Q_cell, A_cell, samples);
    
    [beta_hat, f_min] = Matern_fit(negloglik1, beta_init, lb, ub, @mycon, false);
    rec_beta_hat(rep, :) = beta_hat;
    
end

save(savefile, 'rec_beta_hat');

% run on server
delete(gcp)
