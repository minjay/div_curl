% compute the standard errors by bootstrapping
clear

% run on server
parpool(16)
addpath(genpath('/home/minjay/div_curl'))

load('wind.mat')

savefile = 'boot_NMG.mat';

T = 108;
p = 2;
% number of bootstrap samples
B = 200;

% initial computation
[r, h0_cell] = init_comp_NMG(n, theta, phi, x, y, z);

beta_all = [0.016283 -0.004629 0.002649 0.011662 2.867271 14.523555 0.270156 0.343079 0.720569 0.830202 0.212515 0.194964];
a1 = beta_all(1);
a2 = beta_all(2);
b1 = beta_all(3);
b2 = beta_all(4);
nu = beta_all(5);
a = beta_all(6);
sigma1 = beta_all(7);
sigma2 = beta_all(8);
w1 = beta_all(9);
w2 = beta_all(10);
tau1 = beta_all(11);
tau2 = beta_all(12);
cov_mat_NMG = get_cov_NMG(r, a1, a2, b1, b2, nu, a, h0_cell)+...
    get_cov_Matern_pars(r, sigma1, sigma2, 0, w1, w2, a)+...
    diag(kron(ones(1, n), [tau1^2, tau2^2]));

samples_all = mvnrnd(zeros(p*n, 1), cov_mat, B*T);

samples_all_cell = cell(B, 1);
for rep = 1:B
    samples_all_cell{rep} = samples_all((rep-1)*T+1:rep*T, :);
end

beta_init = beta_all;
% to avoid identifiability problem, set a1>0
lb = [0   -Inf -Inf -Inf 1 0   0   0   0 0 0   0];
ub = [Inf Inf  Inf  Inf  5 Inf Inf Inf 5 5 Inf Inf];

rec_beta_hat = zeros(B, 12);

parfor rep = 1:B
    
    samples = samples_all_cell{rep};
    
    % negative log-likelihood function
    negloglik1 = @(beta_all) negloglik_NMG_Matern_all(beta_all, r, samples, h0_cell);
    
    % fit the model
    [beta_hat, f_min] = Matern_fit(negloglik1, beta_init, lb, ub, [], false);
    rec_beta_hat(rep, :) = beta_hat;
    
end

save(savefile, 'rec_beta_hat');

% run on server
delete(gcp)
