% cokriging
clear

% run on server
parpool(12)
addpath(genpath('/home/minjay/div_curl'))

savefile = 'pred1.mat';

load('wind.mat')

% initial computation
[h_mat, r, P_cell, Q_cell, A_cell] = init_comp(x, y, z, n, theta, phi);

beta_init = [0.029113 0.054503 0.281047 1.757812 2.034312 9.47189 0.210496 0.196141];

lb = [0 0 -1 1 1 0 0 0];
ub = [Inf Inf 1 5 5 Inf Inf Inf];

T = 108;
n_pred = 170;
p = 2;

% get pred and est locations
subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0.1], [0.1 0.01], [0.1 0.01]);
rng('default')
pred_loc = sort(randsample(n, n_pred));
est_loc = setdiff((1:n)', pred_loc);
n_est = n-n_pred;
idx_est = sort([est_loc*p-1; est_loc*p]);
idx_pred = setdiff(1:p*n, idx_est);

lat_est = (pi/2-theta(est_loc))/pi*180;
lon_est = phi(est_loc)/pi*180;
lat_pred = (pi/2-theta(pred_loc))/pi*180;
lon_pred = phi(pred_loc)/pi*180;
h = figure;
subplot(1,1,1)
scatter(lon_est, lat_est, '.');
hold on
scatter(lon_pred, lat_pred, 'r*')
axis equal 
axis tight
xlabel('East longitude')
ylabel('Latitude')

set(h, 'Position', [0, 0, 300, 350]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
[beta_hat, f_min] = Matern_fit(negloglik1, beta_init, lb, ub, @mycon, true);

%beta_hat = [0.028968 0.053642 0.286489 1.802897 2.069642 9.747633 0.212997 0.199817];

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

obs_y = samples(:, idx_pred)';

obs_u = obs_y(1:p:end, :);
obs_v = obs_y(2:p:end, :);
pred_u = pred_y(1:p:end, :);
pred_v = pred_y(2:p:end, :);

save(savefile, 'obs_u', 'obs_v', 'pred_u', 'pred_v');

% plot
h = figure;
subplot(1,2,1)
scatter(obs_u(:), pred_u(:), '.')
hline = refline(1, 0);
set(hline, 'Color', 'g')
set(hline, 'LineWidth', 2)
ls_fit = polyfit(obs_u(:), pred_u(:), 1);
hline = refline(ls_fit);
set(hline, 'Color', 'r')
set(hline, 'LineWidth', 2)
set(hline, 'LineStyle', '--')
axis equal
axis tight
xlabel('Observed')
ylabel('Predicted')
title('U Residual Field')
subplot(1,2,2)
scatter(obs_v(:), pred_v(:), '.')
hline = refline(1, 0);
set(hline, 'Color', 'g')
set(hline, 'LineWidth', 2)
ls_fit = polyfit(obs_v(:), pred_v(:), 1);
hline = refline(ls_fit);
set(hline, 'Color', 'r')
set(hline, 'LineWidth', 2)
set(hline, 'LineStyle', '--')
axis equal
axis tight
xlabel('Observed')
ylabel('Predicted')
title('V Residual Field')
set(h, 'Position', [0, 0, 500, 275]);

% run on server
delete(gcp)
