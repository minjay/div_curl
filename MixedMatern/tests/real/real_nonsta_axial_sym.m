clear

load('wind.mat')

p = 2;

cov_mat_data = cov(samples);
cov_u_data = cov_mat_data(1:p:end, 1:p:end);
cov_v_data = cov_mat_data(2:p:end, 2:p:end);
cov_uv_data = cov_mat_data(1:p:end, 2:p:end);
cov_vu_data = cov_mat_data(2:p:end, 1:p:end);

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
cov_uv = cov_mat(1:p:end, 2:p:end);
cov_vu = cov_mat(2:p:end, 1:p:end);

lat = (pi/2-theta)/pi*180;
lon = phi/pi*180;

% plot covariance of u/v

i1 = 10;
i2 = 10;

[phi_diff, cov_u_data_sub, cov_u_sub, lat1, lat2] = extract_cov(i1, i2, cov_u_data, cov_u, lat, lon);
[phi_diff_sorted, index] = sort(phi_diff);
subplot(2, 2, 1)
plot(phi_diff, cov_u_data_sub, 'bo')
xlabel('\phi_s - \phi_t')
hold on
plot(phi_diff_sorted, cov_u_sub(index), 'r', 'LineWidth', 2)
axis tight
ylim([-0.1 0.3])
title(['Covariance of U on Latitude ', num2str(lat1), '\circ'])
[phi_diff, cov_v_data_sub, cov_v_sub, lat1, lat2] = extract_cov(i1, i2, cov_v_data, cov_v, lat, lon);
[phi_diff_sorted, index] = sort(phi_diff);
subplot(2, 2, 2)
plot(phi_diff, cov_v_data_sub, 'bo')
xlabel('\phi_s - \phi_t')
hold on
plot(phi_diff_sorted, cov_v_sub(index), 'r', 'LineWidth', 2)
axis tight
ylim([-0.1 0.3])
title(['Covariance of V on Latitude ', num2str(lat1), '\circ'])

i1 = 8;
i2 = 10;

[phi_diff, cov_u_data_sub, cov_u_sub, lat1, lat2] = extract_cov(i1, i2, cov_u_data, cov_u, lat, lon);
[phi_diff_sorted, index] = sort(phi_diff);
subplot(2, 2, 3)
plot(phi_diff, cov_u_data_sub, 'bo')
xlabel('\phi_s - \phi_t')
hold on
plot(phi_diff_sorted, cov_u_sub(index), 'r', 'LineWidth', 2)
axis tight
ylim([-0.1 0.3])
title(['Covariance of U on Latitude ', num2str(lat1), '\circ & ', num2str(lat2), '\circ'])
[phi_diff, cov_v_data_sub, cov_v_sub, lat1, lat2] = extract_cov(i1, i2, cov_v_data, cov_v, lat, lon);
[phi_diff_sorted, index] = sort(phi_diff);
subplot(2, 2, 4)
plot(phi_diff, cov_v_data_sub, 'bo')
xlabel('\phi_s - \phi_t')
hold on
plot(phi_diff_sorted, cov_v_sub(index), 'r', 'LineWidth', 2)
axis tight
ylim([-0.1 0.3])
title(['Covariance of V on Latitude ', num2str(lat1), '\circ & ', num2str(lat2), '\circ'])


% plot cross-covariance between u and v

i1 = 10;
i2 = 10;

[phi_diff, cov_uv_data_sub, cov_uv_sub, lat1, lat2] = extract_cov(i1, i2, cov_uv_data, cov_uv, lat, lon);
[phi_diff_sorted, index] = sort(phi_diff);
subplot(2, 2, 1)
plot(phi_diff, cov_uv_data_sub, 'bo')
hold on
plot(phi_diff_sorted, cov_uv_sub(index), 'r', 'LineWidth', 2)
xlabel('\phi_s - \phi_t')
axis tight
ylim([-0.1 0.15])
title(['Cross-covariance of U & V on Latitude ', num2str(lat1), '\circ'])

[phi_diff, cov_vu_data_sub, cov_vu_sub, lat1, lat2] = extract_cov(i1, i2, cov_vu_data, cov_vu, lat, lon);
[phi_diff_sorted, index] = sort(phi_diff);
subplot(2, 2, 2)
plot(phi_diff, cov_vu_data_sub, 'bo')
hold on
plot(phi_diff_sorted, cov_vu_sub(index), 'r', 'LineWidth', 2)
xlabel('\phi_s - \phi_t')
axis tight
ylim([-0.1 0.15])
title(['Cross-covariance of V & U on Latitude ', num2str(lat1), '\circ'])

i1 = 8;
i2 = 10;

[phi_diff, cov_uv_data_sub, cov_uv_sub, lat1, lat2] = extract_cov(i1, i2, cov_uv_data, cov_uv, lat, lon);
[phi_diff_sorted, index] = sort(phi_diff);
subplot(2, 2, 3)
plot(phi_diff, cov_uv_data_sub, 'bo')
hold on
plot(phi_diff_sorted, cov_uv_sub(index), 'r', 'LineWidth', 2)
xlabel('\phi_s - \phi_t')
axis tight
ylim([-0.1 0.15])
title(['Cross-cov of U & V on Lat ', num2str(lat1), '\circ & ', num2str(lat2), '\circ'])

[phi_diff, cov_vu_data_sub, cov_vu_sub, lat1, lat2] = extract_cov(i1, i2, cov_vu_data, cov_vu, lat, lon);
[phi_diff_sorted, index] = sort(phi_diff);
subplot(2, 2, 4)
plot(phi_diff, cov_vu_data_sub, 'bo')
hold on
plot(phi_diff_sorted, cov_vu_sub(index), 'r', 'LineWidth', 2)
xlabel('\phi_s - \phi_t')
axis tight
ylim([-0.1 0.15])
title(['Cross-cov of V & U on Lat ', num2str(lat1), '\circ & ', num2str(lat2), '\circ'])
