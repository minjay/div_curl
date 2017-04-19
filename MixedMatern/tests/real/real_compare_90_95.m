clear
load('wind_90.mat')

p = 2;
lat = (pi/2-theta)/pi*180;
lon = phi/pi*180;

corr_mat_data_90 = corr(samples);
corr_uv_data_90 = corr_mat_data_90(1:p:end, 2:p:end);

load('wind.mat')
corr_mat_data_95 = corr(samples);
corr_uv_data_95 = corr_mat_data_95(1:p:end, 2:p:end);

GREY = [0.6 0.6 0.6];

subplot = @(m,n,p) subtightplot (m, n, p, [0.125 0.075], [0.075 0.05], [0.05 0.01]);

subplot(2, 2, 1)
plot(lat, diag(corr_uv_data_90), 'o', 'Color', GREY, 'MarkerSize', 3)
hold on
y_smooth = smooth(lat, diag(corr_uv_data_90), 0.2, 'loess');
h_emp = plot(lat, y_smooth, 'g', 'LineWidth', 1.5);
axis([min(lat) max(lat) -1 1])
xlabel('Latitude')
title('Cross-correlation of U & V Residual Fields')
legend(h_emp, 'Empirical (90% EOF)')

subplot(2, 2, 2)
plot(lon, diag(corr_uv_data_90), 'o', 'Color', GREY, 'MarkerSize', 3)
hold on
[sorted_lon, index] = sort(lon);
y = diag(corr_uv_data_90);
y_smooth = smooth(sorted_lon, y(index), 0.2, 'loess');
h_emp = plot(sorted_lon, y_smooth, 'g', 'LineWidth', 1.5);
axis([min(lon) max(lon) -1 1])
xlabel('Longitude')
title('Cross-correlation of U & V Residual Fields')
legend(h_emp, 'Empirical (90% EOF)')

subplot(2, 2, 3)
plot(lat, diag(corr_uv_data_95), 'o', 'Color', GREY, 'MarkerSize', 3)
hold on
y_smooth = smooth(lat, diag(corr_uv_data_95), 0.2, 'loess');
h_emp = plot(lat, y_smooth, 'g', 'LineWidth', 1.5);
axis([min(lat) max(lat) -1 1])
xlabel('Latitude')
title('Cross-correlation of U & V Residual Fields')
legend(h_emp, 'Empirical (95% EOF)')

subplot(2, 2, 4)
plot(lon, diag(corr_uv_data_95), 'o', 'Color', GREY, 'MarkerSize', 3)
hold on
[sorted_lon, index] = sort(lon);
y = diag(corr_uv_data_95);
y_smooth = smooth(sorted_lon, y(index), 0.2, 'loess');
h_emp = plot(sorted_lon, y_smooth, 'g', 'LineWidth', 1.5);
axis([min(lon) max(lon) -1 1])
xlabel('Longitude')
title('Cross-correlation of U & V Residual Fields')
legend(h_emp, 'Empirical (95% EOF)')
