clear

load('wind.mat')

t = 1;
p = 2;

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

idx_est = setdiff(1:2*n, [region*2-1; region*2]);

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
Sigma00 = cov_mat_TMM(idx_est, idx_est);
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
obs_u = obs_y(1:p:end);
obs_v = obs_y(2:p:end);
obs_u_region = obs_u(region);
obs_v_region = obs_v(region);

subplot('position', [0.05 0.125 0.3 0.8])
quiver(lon_region, lat_region, obs_u_region, obs_v_region, 'b')
axis equal
axis([min(lon_region)-1 max(lon_region)+1 min(lat_region) max(lat_region)])
rectangle('Position',[87 -25 10 11], 'EdgeColor', 'r', 'LineWidth', 1.5)
title('Observed Wind Field')
xlabel('East longitude')
ylabel('Latitude')

subplot('position', [0.35 0.125 0.3 0.8])
quiver(lon_region, lat_region, pred_u_region_BM, pred_v_region_BM, 'b')
axis equal
axis([min(lon_region)-1 max(lon_region)+1 min(lat_region) max(lat_region)])
rectangle('Position',[87 -25 10 11], 'EdgeColor', 'r', 'LineWidth', 1.5)
title('Predicted Wind Field (PARS-BM)')
xlabel('East longitude')

subplot('position', [0.65 0.125 0.3 0.8])
quiver(lon_region, lat_region, pred_u_region_TMM, pred_v_region_TMM, 'b')
axis equal
axis([min(lon_region)-1 max(lon_region)+1 min(lat_region) max(lat_region)])
rectangle('Position',[87 -25 10 11], 'EdgeColor', 'r', 'LineWidth', 1.5)
title('Predicted Wind Field (TMM)')
xlabel('East longitude')

norm(pred_u_region_BM-obs_u_region)
norm(pred_u_region_TMM-obs_u_region)
norm(pred_v_region_BM-obs_v_region)
norm(pred_v_region_TMM-obs_v_region)