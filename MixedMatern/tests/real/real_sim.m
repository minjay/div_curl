clear

load('wind.mat')

subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0.1], [0.01 0.01], [0.1 0.01]);

p = 2;

% check normality
u1 = samples(1, 1:p:end);
v1 = samples(1, 2:p:end);
h = figure;
subplot(1,3,1)
qqplot(u1)
title('Q-Q Plot of U Residual Field')
axis square
subplot(1,3,2)
qqplot(v1)
title('Q-Q Plot of V Residual Field')
axis square

% Chi-Square Q-Q plot
subplot(1,3,3)
qqchi2([u1' v1'])
xlabel('Chi-Square Quantiles')
ylabel('Squared Mahalanobis Distances')
title('Chi-Square Q-Q Plot')
axis square

set(h, 'Position', [0, 0, 600, 200]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cov_mat_data = cov(samples);
cov_u_data = cov_mat_data(1:p:end, 1:p:end);
cov_v_data = cov_mat_data(2:p:end, 2:p:end);

% initial computation
[h_mat, r, P_cell, Q_cell, A_cell] = init_comp(x, y, z, n, theta, phi);
[r, h0_cell] = init_comp_NMG(n, theta, phi, x, y, z);

idx_upper = find(triu(ones(n)));
r_upper = r(idx_upper);
r_upper = 2*asin(r_upper/2)*6371;
[r_upper_sort, r_upper_index] = sort(r_upper);
cov_u_data_upper = cov_u_data(idx_upper);
cov_v_data_upper = cov_v_data(idx_upper);

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

cov_u_upper = cov_u(idx_upper);
cov_v_upper = cov_v(idx_upper);
% number of bins
nb = 100;

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
cov_uv_BM = cov_mat_BM(1:p:end, 2:p:end);
cov_u_upper_BM = cov_u_BM(idx_upper);
cov_v_upper_BM = cov_v_BM(idx_upper);

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
cov_uv_NMG = cov_mat_NMG(1:p:end, 2:p:end);
cov_vu_NMG = cov_mat_NMG(2:p:end, 1:p:end);
cov_u_upper_NMG = cov_u_NMG(idx_upper);
cov_v_upper_NMG = cov_v_NMG(idx_upper);

width = 75;
h = figure;

subplot = @(m,n,p) subtightplot (m, n, p, [0.125 0.075], [0.075 0.05], [0.05 0.01]);

GREY = [0.6 0.6 0.6];

subplot(2,2,1)
[X_MED,Y_MED,~,~] = binned_plot(r_upper, cov_u_data_upper);
r_max_index = find(r_upper_sort<=X_MED(end), 1, 'last');
plot(X_MED, Y_MED, '.', 'Color', GREY, 'MarkerSize', 8)
lh = boxplot_curve(r_upper, cov_u_upper, nb, 'b', 0.5);
ph = plot(r_upper_sort(1:r_max_index), cov_u_upper_BM(r_upper_index(1:r_max_index)), 'k-.', 'LineWidth', 1.5);
lh2 = boxplot_curve2(r_upper, cov_u_upper_NMG, nb, 'r', 0.5, width);
axis tight
xlabel('Great-circle Distance (km)')
title('Covariance of U Residual Field')
legend([lh ph lh2], {'TMM', 'PARS-BM', 'NBG'})

subplot(2,2,2)
[X_MED,Y_MED,~,~] = binned_plot(r_upper, cov_v_data_upper);
plot(X_MED, Y_MED, '.', 'Color', GREY, 'MarkerSize', 8)
lh = boxplot_curve(r_upper, cov_v_upper, nb, 'b', 0.5);
ph = plot(r_upper_sort(1:r_max_index), cov_v_upper_BM(r_upper_index(1:r_max_index)), 'k-.', 'LineWidth', 1.5);
lh2 = boxplot_curve2(r_upper, cov_v_upper_NMG, nb, 'r', 0.5, width);
axis tight
xlabel('Great-circle Distance (km)')
title('Covariance of V Residual Field')
legend([lh ph lh2], {'TMM', 'PARS-BM' 'NBG'})

lat = (pi/2-theta)/pi*180;
lon = phi/pi*180;

% compute the correlations
corr_mat_data = corr(samples);
corr_uv_data = corr_mat_data(1:p:end, 2:p:end);

corr_uv_BM = diag(1./sqrt(diag(cov_u_BM)))*cov_uv_BM*diag(1./sqrt(diag(cov_v_BM)));
corr_uv_NMG = diag(1./sqrt(diag(cov_u_NMG)))*cov_uv_NMG*diag(1./sqrt(diag(cov_v_NMG)));

subplot(2,2,3)
plot(lat, diag(corr_uv_data), '.', 'Color', GREY, 'MarkerSize', 8)
hold on
y_smooth = smooth(lat, diag(corr_uv_data), 0.2, 'loess');
h_emp = plot(lat, y_smooth, 'g', 'LineWidth', 1.5);
axis([min(lat) max(lat) -1 1])
hline = refline(0, 0);
set(hline,'Color','b')
set(hline,'LineWidth', 1.5)
ph = plot(lat, diag(corr_uv_BM), 'k-.', 'LineWidth', 1.5);
ph2 = plot(lat, diag(corr_uv_NMG), 'r--', 'LineWidth', 1.5);
xlabel('Latitude (degree)')
title('Cross-correlation of U & V Residual Fields')
legend([h_emp hline ph ph2], {'Empirical','TMM', 'PARS-BM', 'NBG'})

subplot(2,2,4)
plot(lon, diag(corr_uv_data), '.', 'Color', GREY, 'MarkerSize', 8)
hold on
[sorted_lon, index] = sort(lon);
y = diag(corr_uv_data);
y_smooth = smooth(sorted_lon, y(index), 0.2, 'loess');
h_emp = plot(sorted_lon, y_smooth, 'g', 'LineWidth', 1.5);
axis([min(lon) max(lon) -1 1])
hline = refline(0, 0);
set(hline,'Color','b')
set(hline,'LineWidth', 1.5)
ph = plot(sorted_lon, diag(corr_uv_BM(index, index)), 'k-.', 'LineWidth', 1.5);
ph2 = plot(sorted_lon, diag(corr_uv_NMG(index, index)), 'r.');
xlabel('Longitude (degree)')
title('Cross-correlation of U & V Residual Fields')
legend([h_emp hline ph ph2], {'Empirical','TMM', 'PARS-BM', 'NBG'})
