clear

load('wind.mat')

p = 2;

cov_mat_data = cov(samples);
cov_u_data = cov_mat_data(1:p:end, 1:p:end);
cov_v_data = cov_mat_data(2:p:end, 2:p:end);

% initial computation
[h_mat, r, P_cell, Q_cell, A_cell] = init_comp(x, y, z, n, theta, phi);
[r, h0_cell] = init_comp_NMG(n, theta, phi, x, y, z);

% specify parameters
beta_all = [0.029113 0.054503 0.281047 1.757812 2.034312 9.47189 0.210496 0.196141];
beta = beta_all(1:6);
sigma1 = beta(1);
sigma2 = beta(2);
nu1 = beta(4);
nu2 = beta(5);
a = beta(6);
tau1 = beta_all(7);
tau2 = beta_all(8);

% compute the signal-to-noise ratios
signal = sigma1^2*a^2/2/(nu1-1)+sigma2^2*a^2/2/(nu2-1);
prop1 = signal/tau1^2;
prop2 = signal/tau2^2;

% get coef and bessel
[coef, bessel] = get_coef_bessel(beta, r);

% get cov mat
cov_mat = get_cov(h_mat, r, P_cell, Q_cell, A_cell, @Matern_mix,...
    beta, coef, bessel)+diag(kron(ones(1, n), [tau1^2, tau2^2]));

cov_u = cov_mat(1:p:end, 1:p:end);
cov_v = cov_mat(2:p:end, 2:p:end);

% get cov mat for BM
param_BM = [0.21810846 0.203180249 0.058196367 1.23929461 1.1322247 0.184217146 0.156812964 -0.07977680];
beta_all(1:2) = param_BM(6:7);
beta_all(3) = param_BM(8);
beta_all(4:5) = param_BM(4:5);
beta_all(6) = param_BM(3);
beta_all(7:8) = param_BM(1:2);
sigma1 = sqrt(beta_all(1));
sigma2 = sqrt(beta_all(2));
rho12 = beta_all(3);
nu1 = beta_all(4);
nu2 = beta_all(5);
a = 1/beta_all(6);
tau1 = beta_all(7);
tau2 = beta_all(8);
cov_mat_BM = get_cov_Matern_pars(r, sigma1, sigma2, rho12, nu1, nu2, a)+...
    diag(kron(ones(1, n), [tau1^2, tau2^2]));
cov_u_BM = cov_mat_BM(1:p:end, 1:p:end);
cov_v_BM = cov_mat_BM(2:p:end, 2:p:end);

% get cov mat for NMG
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
cov_u_NMG = cov_mat_NMG(1:p:end, 1:p:end);
cov_v_NMG = cov_mat_NMG(2:p:end, 2:p:end);

lat = (pi/2-theta)/pi*180;
lon = phi/pi*180;
[lon_sorted, index] = sort(lon);

subplot = @(m,n,p) subtightplot (m, n, p, [0.125 0.075], [0.075 0.05], [0.05 0.02]);
GREY = [0.6 0.6 0.6];

% compare true and fitted variance at each location
% plot against latitudes and longitudes 
subplot(2, 2, 1)
plot(lat, diag(cov_u_data), '.', 'Color', GREY, 'MarkerSize', 8)
hold on
y_smooth = smooth(lat, diag(cov_u_data), 0.2, 'loess');
h_emp = plot(lat, y_smooth, 'g', 'LineWidth', 1.5);
ph1 = plot(lat, diag(cov_u), 'b', 'LineWidth', 1.5);
ph2 = plot(lat, diag(cov_u_BM), 'k-.', 'LineWidth', 1.5);
ph3 = plot(lat, diag(cov_u_NMG), 'r--', 'LineWidth', 1.5);
axis tight
title('Variance of U Residual Field')
ylim([0 1])
xlabel('Latitude')
legend([h_emp ph1 ph2 ph3], {'Empirical', 'TMM', 'PARS-BM', 'NBG'})

subplot(2, 2, 2)
plot(lon, diag(cov_u_data), '.', 'Color', GREY, 'MarkerSize', 8)
hold on
y_smooth = smooth(lon_sorted, diag(cov_u_data(index, index)), 0.2, 'loess');
h_emp = plot(lon_sorted, y_smooth, 'g', 'LineWidth', 1.5);
ph1 = plot(lon_sorted, diag(cov_u(index, index)), 'b', 'LineWidth', 1.5);
ph2 = plot(lon_sorted, diag(cov_u_BM(index, index)), 'k-.', 'LineWidth', 1.5);
ph3 = plot(lon_sorted, diag(cov_u_NMG(index, index)), 'r.');
axis tight
title('Variance of U Residual Field')
ylim([0 1])
xlabel('Longitude')
legend([h_emp ph1 ph2 ph3], {'Empirical', 'TMM', 'PARS-BM', 'NBG'})

subplot(2, 2, 3)
plot(lat, diag(cov_v_data), '.', 'Color', GREY, 'MarkerSize', 8)
hold on
y_smooth = smooth(lat, diag(cov_v_data), 0.2, 'loess');
h_emp = plot(lat, y_smooth, 'g', 'LineWidth', 1.5);
ph1 = plot(lat, diag(cov_v), 'b', 'LineWidth', 1.5);
ph2 = plot(lat, diag(cov_v_BM), 'k-.', 'LineWidth', 1.5);
ph3 = plot(lat, diag(cov_v_NMG), 'r--', 'LineWidth', 1.5);
axis tight
title('Variance of V Residual Field')
ylim([0 1])
xlabel('Latitude')
legend([h_emp ph1 ph2 ph3], {'Empirical', 'TMM', 'PARS-BM', 'NBG'})

subplot(2, 2, 4)
plot(lon, diag(cov_v_data), '.', 'Color', GREY, 'MarkerSize', 8)
hold on
y_smooth = smooth(lon_sorted, diag(cov_v_data(index, index)), 0.2, 'loess');
h_emp = plot(lon_sorted, y_smooth, 'g', 'LineWidth', 1.5);
ph1 = plot(lon_sorted, diag(cov_v(index, index)), 'b', 'LineWidth', 1.5);
ph2 = plot(lon_sorted, diag(cov_v_BM(index, index)), 'k-.', 'LineWidth', 1.5);
ph3 = plot(lon_sorted, diag(cov_v_NMG(index, index)), 'r.');
axis tight
title('Variance of V Residual Field')
ylim([0 1])
xlabel('Longitude')
legend([h_emp ph1 ph2 ph3], {'Empirical', 'TMM', 'PARS-BM', 'NBG'})
