% compute the standard errors by bootstrapping
% on Hannan
clear

% run on server
parpool(16)
addpath(genpath('/home/minjay/div_curl'))

load('wind.mat')
load('samples_NMG.mat')

savefile = 'boot_NMG2.mat';

T = 108;
p = 2;
% number of bootstrap samples
B = 200;

% initial computation
[r, h0_cell] = init_comp_NMG(n, theta, phi, x, y, z);

beta_all = [0.016283 -0.004629 0.002649 0.011662 2.867271 14.523555 0.270156 0.343079 0.720569 0.830202 0.212515 0.194964];

beta_init = beta_all;
% to avoid identifiability problem, set a1>0
lb = [0   -Inf -Inf -Inf 1 0   0   0   0 0 0   0];
ub = [Inf Inf  Inf  Inf  5 Inf Inf Inf 5 5 Inf Inf];

rec_beta_hat = zeros(B, 12);

parfor rep = 51:100
    
    samples = samples_all_cell{rep};
    
    % negative log-likelihood function
    negloglik1 = @(beta_all) negloglik_NMG_Matern_all(beta_all, r, samples, h0_cell);
    
    % fit the model
    [beta_hat, f_min] = Matern_fit_Nelder(negloglik1, beta_init);
    rec_beta_hat(rep, :) = beta_hat;
    
end

save(savefile, 'rec_beta_hat');

% run on server
delete(gcp)
