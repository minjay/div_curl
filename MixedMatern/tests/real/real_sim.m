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

idx_upper = find(triu(ones(n)));
r_upper = r(idx_upper);
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

h = figure;

subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0.05], [0.05 0.05], [0.05 0.01]);

r_upper = 2*asin(r_upper/2)*6371;

subplot(2,2,1)
[X_MED,Y_MED,~,~] = binned_plot(r_upper, cov_u_data_upper);
scatter(X_MED, Y_MED)
boxplot_curve(r_upper, cov_u_upper, nb, 'r')
axis square
axis tight
xlabel('Great-circle Distance (km)')
title('Covariance of U Residual Field')
subplot(2,2,2)
[X_MED,Y_MED,~,~] = binned_plot(r_upper, cov_v_data_upper);
scatter(X_MED, Y_MED)
boxplot_curve(r_upper, cov_v_upper, nb, 'r')
axis square
axis tight
xlabel('Great-circle Distance (km)')
title('Covariance of V Residual Field')

lat = (pi/2-theta)/pi*180;
lon = phi/pi*180;

% compute the correlation
corr_mat_data = corr(samples);
corr_uv_data = corr_mat_data(1:p:end, 2:p:end);

subplot(2,2,3)
scatter(lat, diag(corr_uv_data))
hold on
[X_MED,Y_MED,~,~] = binned_plot(lat, diag(corr_uv_data));
h_emp = plot(X_MED, Y_MED, 'g', 'LineWidth', 2);
axis square
axis([min(lat) max(lat) -1 1])
hline = refline(0, 0);
set(hline,'Color','r')
set(hline,'LineWidth', 2)
set(hline,'LineStyle', '--')
xlabel('Latitude')
title('Cross-correlation of U & V Residual Fields')
legend([h_emp hline], {'Empirical','Fitted'})

subplot(2,2,4)
scatter(lon, diag(corr_uv_data))
hold on
[X_MED,Y_MED,~,~] = binned_plot(lon, diag(corr_uv_data));
h_emp = plot(X_MED, Y_MED, 'g', 'LineWidth', 2);
axis square
axis([min(lon) max(lon) -1 1])
hline = refline(0, 0);
set(hline,'Color','r')
set(hline,'LineWidth', 2)
set(hline,'LineStyle', '--')
xlabel('Longitude')
title('Cross-correlation of U & V Residual Fields')
legend([h_emp hline], {'Empirical','Fitted'})

set(h, 'Position', [0, 0, 600, 600]);
