clear

load('wind.mat')

t = 1;
p = 2;

% initial computation
[h_mat, r, P_cell, Q_cell, A_cell] = init_comp(x, y, z, n, theta, phi);
[r, h0_cell] = init_comp_NMG(n, theta, phi, x, y, z);

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

% NMG
param_NMG = [0.016283 -0.004629 0.002649 0.011662 2.867271 14.523555 0.270156 0.343079 0.720569 0.830202 0.212515 0.194964];
a1 = param_NMG(1);
a2 = param_NMG(2);
b1 = param_NMG(3);
b2 = param_NMG(4);
nu = param_NMG(5);
a = param_NMG(6);
sigma1 = param_NMG(7);
sigma2 = param_NMG(8);
w1 = param_NMG(9);
w2 = param_NMG(10);
tau1 = param_NMG(11);
tau2 = param_NMG(12);
cov_mat_NMG = get_cov_NMG(r, a1, a2, b1, b2, nu, a, h0_cell)+...
    get_cov_Matern_pars(r, sigma1, sigma2, 0, w1, w2, a)+...
    diag(kron(ones(1, n), [tau1^2, tau2^2]));

% follow the cokriging formula
Sigma00 = cov_mat_NMG(idx_est, idx_est);
est_y = samples(t, idx_est)';
tmp = Sigma00\est_y;

SigmaP0 = cov_mat_NMG(:, idx_est);
pred_y = SigmaP0*tmp;

pred_u = pred_y(1:p:end);
pred_v = pred_y(2:p:end);

pred_u_region_NMG = pred_u(region);
pred_v_region_NMG = pred_v(region);

% TMM
beta_hat = [0.029113 0.054503 0.281047 1.757812 2.034312 9.47189 0.210496 0.196141];

beta = beta_hat(1:6);
tau1 = beta_hat(7);
tau2 = beta_hat(8);

[coef, bessel] = get_coef_bessel(beta, r);

% get cov mat
cov_mat_TMM = get_cov(h_mat, r, P_cell, Q_cell, A_cell, @Matern_mix,...
    beta, coef, bessel)+diag(kron(ones(1, n), [tau1^2, tau2^2]));

% get Cov(curl-free(s), y(t)) and Cov(div-free(s), y(t))
cov_mat_curl = get_cov_full(h_mat, r, P_cell, Q_cell, A_cell, @Matern_curl,...
    beta, coef, bessel);
cov_mat_div = get_cov_full(h_mat, r, P_cell, Q_cell, A_cell, @Matern_div,...
    beta, coef, bessel);

% cov_mat_TMM should be equal to cov_mat_curl+cov_mat_div
tmp = cov_mat_curl+cov_mat_div+diag(kron(ones(1, n), [tau1^2, tau2^2]));
max(max((abs(cov_mat_TMM-tmp))))

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

% predict curl-free component    
SigmaP0_curl = cov_mat_curl(:, idx_est);
pred_y_curl = SigmaP0_curl*tmp;

pred_u_curl = pred_y_curl(1:p:end);
pred_v_curl = pred_y_curl(2:p:end);
pred_u_region_curl = pred_u_curl(region);
pred_v_region_curl = pred_v_curl(region);

% predict div-free component
SigmaP0_div = cov_mat_div(:, idx_est);
pred_y_div = SigmaP0_div*tmp;

pred_u_div = pred_y_div(1:p:end);
pred_v_div = pred_y_div(2:p:end);
pred_u_region_div = pred_u_div(region);
pred_v_region_div = pred_v_div(region);

% sanity checking
pred_u_region_curl_div = pred_u_region_curl+pred_u_region_div;
max(max(abs(pred_u_region_curl_div-pred_u_region_TMM)))
pred_v_region_curl_div = pred_v_region_curl+pred_v_region_div;
max(max(abs(pred_v_region_curl_div-pred_v_region_TMM)))

% plot
obs_y = samples(t, :)';
obs_u = obs_y(1:p:end);
obs_v = obs_y(2:p:end);
obs_u_region = obs_u(region);
obs_v_region = obs_v(region);

subplot = @(m,n,p) subtightplot (m, n, p, [0.125 0.05], [0.075 0.05], [0.05 0.02]);

subplot(2, 3, 1)
% note that "quiver" automatically scales the arrows to fit within the grid
quiver(lon_region, lat_region, obs_u_region, obs_v_region, 'b')
axis equal
axis([min(lon_region)-1 max(lon_region)+1 min(lat_region) max(lat_region)])
rectangle('Position',[87 -25 10 11], 'EdgeColor', 'r', 'LineWidth', 1.5)
title('Observed')
xlabel('East longitude')
ylabel('Latitude')

subplot(2, 3, 2)
quiver(lon_region, lat_region, pred_u_region_BM, pred_v_region_BM, 'b')
axis equal
axis([min(lon_region)-1 max(lon_region)+1 min(lat_region) max(lat_region)])
rectangle('Position',[87 -25 10 11], 'EdgeColor', 'r', 'LineWidth', 1.5)
title('Predicted (PARS-BM)')
xlabel('East longitude')

subplot(2, 3, 3)
quiver(lon_region, lat_region, pred_u_region_NMG, pred_v_region_NMG, 'b')
axis equal
axis([min(lon_region)-1 max(lon_region)+1 min(lat_region) max(lat_region)])
rectangle('Position',[87 -25 10 11], 'EdgeColor', 'r', 'LineWidth', 1.5)
title('Predicted (NBG)')
xlabel('East longitude')

norm(pred_u_region_BM-obs_u_region)
norm(pred_u_region_TMM-obs_u_region)
norm(pred_v_region_BM-obs_v_region)
norm(pred_v_region_TMM-obs_v_region)

subplot(2, 3, 4)
quiver(lon_region, lat_region, pred_u_region_TMM, pred_v_region_TMM, 'b')
axis equal
axis([min(lon_region)-1 max(lon_region)+1 min(lat_region) max(lat_region)])
rectangle('Position',[87 -25 10 11], 'EdgeColor', 'r', 'LineWidth', 1.5)
title('Predicted (TMM)')
xlabel('East longitude')
ylabel('Latitude')

subplot(2, 3, 5)
quiver(lon_region, lat_region, pred_u_region_curl, pred_v_region_curl, 'b')
axis equal
axis([min(lon_region)-1 max(lon_region)+1 min(lat_region) max(lat_region)])
title('Predicted Curl-free (TMM)')
xlabel('East longitude')

subplot(2, 3, 6)
quiver(lon_region, lat_region, pred_u_region_div, pred_v_region_div, 'b')
axis equal
axis([min(lon_region)-1 max(lon_region)+1 min(lat_region) max(lat_region)])
title('Predicted Div-free (TMM)')
xlabel('East longitude')
