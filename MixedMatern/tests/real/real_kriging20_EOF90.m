% cokriging 20 times
% work on residuals after removing 90% EOFs
clear

% run on server
parpool(12)
addpath(genpath('/home/minjay/div_curl'))

savefile = 'pred_err_90.mat';

load('wind_90.mat')

% initial computation
[h_mat, r, P_cell, Q_cell, A_cell] = init_comp(x, y, z, n, theta, phi);

beta_init = [0.029113 0.054503 0.281047 1.757812 2.034312 9.47189 0.210496 0.196141];

lb = [0 0 -1 1 1 0 0 0];
ub = [Inf Inf 1 5 5 Inf Inf Inf];

T = 108;
n_pred = n/2;
p = 2;
B = 20;

MSPE_u = zeros(B, 1);
MSPE_v = zeros(B, 1);
MAE_u = zeros(B, 1);
MAE_v = zeros(B, 1);

LogS_u = zeros(B, 1);
LogS_v = zeros(B, 1);
CRPS_u = zeros(B, 1);
CRPS_v = zeros(B, 1);

rng('default')

n_est = n-n_pred;
rec_pred_loc = zeros(B, n_pred);
rec_est_loc = zeros(B, n_est);
rec_idx_est = zeros(B, p*n_est);
rec_idx_pred = zeros(B, p*n_pred);

theta_m = (min(theta)+max(theta))/2;
phi_m = (min(phi)+max(phi))/2;
width = 10/180*pi;
theta_l = theta_m-2*width;
theta_r = theta_m+2*width;
phi_l = phi_m-width;
phi_r = phi_m+width;

region = find(theta>=theta_l & theta<=theta_r & phi>=phi_l & phi<=phi_r);
pop = setdiff(1:n, region);

for i = 1:B
    rec_est_loc(i, :) = sort(randsample(pop, n_est));
    rec_pred_loc(i, :) = setdiff(1:n, rec_est_loc(i, :));

    rec_idx_est(i, :) = sort([rec_est_loc(i, :)*p-1 rec_est_loc(i, :)*p]);
    rec_idx_pred(i, :) = setdiff(1:p*n, rec_idx_est(i, :));
end

save('pred_loc_90.mat', 'rec_pred_loc', 'rec_est_loc', 'rec_idx_est', 'rec_idx_pred')

parfor rep = 1:B
    
    % get pred and est locations
    pred_loc = rec_pred_loc(rep, :);
    est_loc = rec_est_loc(rep, :);

    idx_est = rec_idx_est(rep, :);
    idx_pred = rec_idx_pred(rep, :);

    % compute the subsets
    h_mat_sub = h_mat(est_loc, est_loc);
    r_sub = r(est_loc, est_loc);
    P_cell_sub = P_cell(est_loc);
    Q_cell_sub = Q_cell(est_loc);
    A_cell_sub = A_cell(est_loc);
    samples_sub = samples(:, idx_est);

    % negative log-likelihood function
    negloglik1 = @(beta_all) negloglik(beta_all, h_mat_sub, r_sub, P_cell_sub, Q_cell_sub, A_cell_sub, samples_sub);

    % fit the model
    [beta_hat, f_min] = Matern_fit(negloglik1, beta_init, lb, ub, @mycon, false);

    beta = beta_hat(1:6);
    tau1 = beta_hat(7);
    tau2 = beta_hat(8);

    [coef, bessel] = get_coef_bessel(beta, r);

    % get cov mat
    cov_mat = get_cov(h_mat, r, P_cell, Q_cell, A_cell, @Matern_mix,...
        beta, coef, bessel)+diag(kron(ones(1, n), [tau1^2, tau2^2]));

    % follow the cokriging formula
    Sigma00 = cov_mat(idx_est, idx_est);
    est_y = samples(:, idx_est)';
    tmp = Sigma00\est_y;
    
    SigmaP0 = cov_mat(idx_pred, idx_est);
    pred_y = SigmaP0*tmp;
    
    SigmaPP = cov_mat(idx_pred, idx_pred);
    Sigma0P = SigmaP0';
    var_pred_y = diag(SigmaPP-SigmaP0*(Sigma00\Sigma0P));
    var_pred_y_u = var_pred_y(1:p:end);
    var_pred_y_v = var_pred_y(2:p:end);

    obs_y = samples(:, idx_pred)';

    obs_u = obs_y(1:p:end, :);
    obs_v = obs_y(2:p:end, :);
    pred_u = pred_y(1:p:end, :);
    pred_v = pred_y(2:p:end, :);
    
    diff_u = obs_u-pred_u;
    diff_u = diff_u(:);
    diff_v = obs_v-pred_v;
    diff_v = diff_v(:);
    
    negloglik_u = zeros(T, 1);
    negloglik_v = zeros(T, 1);
    CRPS_t_u = zeros(T, 1);
    CRPS_t_v = zeros(T, 1);
    for t = 1:T
        negloglik_u(t) = mean(-log(normpdf(obs_u(:, t), pred_u(:, t), sqrt(var_pred_y_u))));
        negloglik_v(t) = mean(-log(normpdf(obs_v(:, t), pred_v(:, t), sqrt(var_pred_y_v))));
        CRPS_t_u(t) = mean(CRPS(obs_u(:, t), pred_u(:, t), var_pred_y_u));
        CRPS_t_v(t) = mean(CRPS(obs_v(:, t), pred_v(:, t), var_pred_y_v));
    end
    
    MSPE_u(rep) = mean(diff_u.^2);
    MSPE_v(rep) = mean(diff_v.^2);
    MAE_u(rep) = mean(abs(diff_u));  
    MAE_v(rep) = mean(abs(diff_v));
    LogS_u(rep) = mean(negloglik_u);
    LogS_v(rep) = mean(negloglik_v);
    CRPS_u(rep) = mean(CRPS_t_u);
    CRPS_v(rep) = mean(CRPS_t_v);
    
end

save(savefile, 'MSPE_u', 'MSPE_v', 'MAE_u', 'MAE_v', 'LogS_u', 'LogS_v', 'CRPS_u', 'CRPS_v');

% run on server
delete(gcp)
