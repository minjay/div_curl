clear

load('wind.mat')
load('pred_loc.mat')

rep = 1;
t = 1;
p = 2;

est_loc = rec_est_loc(rep, :);
pred_loc = rec_pred_loc(rep, :);
idx_est = rec_idx_est(rep, :);
idx_pred = rec_idx_pred(rep, :);

% initial computation
[h_mat, r, P_cell, Q_cell, A_cell] = init_comp(x, y, z, n, theta, phi);

% get region
theta_m = (min(theta)+max(theta))/2;
phi_m = (min(phi)+max(phi))/2;
width = 10/180*pi;
theta_l = theta_m-2*width;
theta_r = theta_m+2*width;
phi_l = phi_m-width;
phi_r = phi_m+width;

region = find(theta>=theta_l & theta<=theta_r & phi>=phi_l & phi<=phi_r);

% convert to lat and lon
lat = (pi/2-theta)/pi*180;
lon = phi/pi*180;
lat_region = lat(region);
lon_region = lon(region);

% PARS-BM
sigma1 = sqrt(0.184);
sigma2 = sqrt(0.157);
rho12 = -0.08;
nu1 = 1.239;
nu2 = 1.132;
a = 1/0.058;
tau1 = 0.218;
tau2 = 0.203;
cov_mat_BM = get_cov_Matern_pars(r, sigma1, sigma2, rho12, nu1, nu2, a)+...
    diag(kron(ones(1, n), [tau1^2, tau2^2]));

% follow the cokriging formula
Sigma00 = cov_mat_BM(idx_est, idx_est);
est_y = samples(t, idx_est)';
tmp = Sigma00\est_y;

SigmaP0 = cov_mat_BM(:, idx_est);
pred_y = SigmaP0*tmp;

pred_u = pred_y(1:p:end);
pred_v = pred_y(2:p:end);

pred_u_region_BM = pred_u(region);
pred_v_region_BM = pred_v(region);

% TMM
beta_hat = [0.029113 0.054503 0.281047 1.757812 2.034312 9.47189 0.210496 0.196141];

beta = beta_hat(1:6);
tau1 = beta_hat(7);
tau2 = beta_hat(8);

[coef, bessel] = get_coef_bessel(beta, r);

% get cov mat
cov_mat_TMM = get_cov(h_mat, r, P_cell, Q_cell, A_cell, @Matern_mix,...
    beta, coef, bessel)+diag(kron(ones(1, n), [tau1^2, tau2^2]));

% follow the cokriging formula
Sigma00 = cov_mat_TMM(idx_est, idx_est)+diag(kron(ones(1, n/2), [tau1^2, tau2^2]));
est_y = samples(t, idx_est)';
tmp = Sigma00\est_y;
    
SigmaP0 = cov_mat_TMM(:, idx_est);
pred_y = SigmaP0*tmp;

pred_u = pred_y(1:p:end);
pred_v = pred_y(2:p:end);
pred_u_region_TMM = pred_u(region);
pred_v_region_TMM = pred_v(region);

% plot
obs_y = samples(t, :)';
obs_u = obs_y(1:p:end, :);
obs_v = obs_y(2:p:end, :);
obs_u_region = obs_u(region);
obs_v_region = obs_v(region);

subplot(1, 3, 1)
quiver(lon_region, lat_region, obs_u_region, obs_v_region)
axis equal
axis tight
rectangle('Position',[87 -25 12 11], 'EdgeColor','r')
title('Observed Winds')

subplot(1, 3, 2)
quiver(lon_region, lat_region, pred_u_region_BM, pred_v_region_BM)
axis equal
axis tight
rectangle('Position',[87 -25 12 11], 'EdgeColor','r')
title('Predicted Winds (PARS-BM)')

subplot(1, 3, 3)
quiver(lon_region, lat_region, pred_u_region_TMM, pred_v_region_TMM)
axis equal
axis tight
rectangle('Position',[87 -25 12 11], 'EdgeColor','r')
title('Predicted Winds (TMM)')