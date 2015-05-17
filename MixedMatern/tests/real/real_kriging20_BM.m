% cokriging 20 times
% parsimonious bivariate Matern model
clear

p = 2;
B = 20;

load('wind.mat')
load('pred_loc.mat')

% parameters are estimated by R package RandomFields
load('param_kriging.mat')

% initial computation
[h_mat, r, P_cell, Q_cell, A_cell] = init_comp(x, y, z, n, theta, phi);

MSPE_u = zeros(B, 1);
MSPE_v = zeros(B, 1);
MAE_u = zeros(B, 1);
MAE_v = zeros(B, 1);

beta_all = zeros(1, 8);

for rep = 1:B
    
    idx_est = rec_idx_est(rep, :);
    idx_pred = rec_idx_pred(rep, :);
    
    beta_all(1:2) = param_BM(rep, 6:7);
    beta_all(3) = param_BM(rep, 8);
    beta_all(4:5) = param_BM(rep, 4:5);
    beta_all(6) = param_BM(rep, 3);
    beta_all(7:8) = param_BM(rep, 1:2);
    sigma1 = sqrt(beta_all(1));
    sigma2 = sqrt(beta_all(2));
    rho12 = beta_all(3);
    nu1 = beta_all(4);
    nu2 = beta_all(5);
    a = 1/beta_all(6);
    tau1 = beta_all(7);
    tau2 = beta_all(8);
    cov_mat = get_cov_Matern_pars(r, sigma1, sigma2, rho12, nu1, nu2, a)+...
        diag(kron(ones(1, n), [tau1^2, tau2^2]));
    
    % follow the cokriging formula
    Sigma00 = cov_mat(idx_est, idx_est);
    est_y = samples(:, idx_est)';
    tmp = Sigma00\est_y;
    
    SigmaP0 = cov_mat(idx_pred, idx_est);
    pred_y = SigmaP0*tmp;

    obs_y = samples(:, idx_pred)';

    obs_u = obs_y(1:p:end, :);
    obs_v = obs_y(2:p:end, :);
    pred_u = pred_y(1:p:end, :);
    pred_v = pred_y(2:p:end, :);
    
    diff_u = obs_u-pred_u;
    diff_u = diff_u(:);
    diff_v = obs_v-pred_v;
    diff_v = diff_v(:);
    
    MSPE_u(rep) = mean(diff_u.^2);
    MSPE_v(rep) = mean(diff_v.^2);
    MAE_u(rep) = mean(abs(diff_u));  
    MAE_v(rep) = mean(abs(diff_v));
    
end

% MSPE_u
median(MSPE_u)
max(MSPE_u)
min(MSPE_u)

% MSPE_v
median(MSPE_v)
max(MSPE_v)
min(MSPE_v)

% MAE_u
median(MAE_u)
max(MAE_u)
min(MAE_u)

% MAE_v
median(MAE_v)
max(MAE_v)
min(MAE_v)
