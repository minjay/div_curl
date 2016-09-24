clear

load('wind.mat')

p = 2;

cov_mat_data = cov(samples);
cov_u_data = cov_mat_data(1:p:end, 1:p:end);
cov_v_data = cov_mat_data(2:p:end, 2:p:end);

% initial computation
[h_mat, r, P_cell, Q_cell, A_cell] = init_comp(x, y, z, n, theta, phi);

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

lat = (pi/2-theta)/pi*180;
lon = phi/pi*180;

% compare true and fitted variance at each location
% plot against latitudes and longitudes 
subplot(2, 2, 1)
plot(lat, diag(cov_u_data), 'bo')
hold on
plot(lat, diag(cov_u), 'r', 'LineWidth', 2)
axis tight
title('Variance of U Residual Field')
ylim([0 1])
xlabel('Latitude')
subplot(2, 2, 2)
plot(lon, diag(cov_u_data), 'bo')
hold on
plot(lon, diag(cov_u), 'r', 'LineWidth', 2)
axis tight
title('Variance of U Residual Field')
ylim([0 1])
xlabel('Longitude')
subplot(2, 2, 3)
plot(lat, diag(cov_v_data), 'bo')
hold on
plot(lat, diag(cov_v), 'r', 'LineWidth', 2)
axis tight
title('Variance of V Residual Field')
ylim([0 1])
xlabel('Latitude')
subplot(2, 2, 4)
plot(lon, diag(cov_v_data), 'bo')
hold on
plot(lon, diag(cov_v), 'r', 'LineWidth', 2)
axis tight
title('Variance of V Residual Field')
ylim([0 1])
xlabel('Longitude')