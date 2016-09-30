clear

% read data
uwnd_all = ncread('uas_QuikSCAT_L2B_v20110531_199908-200910.nc', 'uas');
% lat and lon (in radian)
% lat in [-pi/2, pi/2]
% lon in [0, 2*pi]
lat = ncread('uas_QuikSCAT_L2B_v20110531_199908-200910.nc', 'lat')/180*pi;
lon = ncread('uas_QuikSCAT_L2B_v20110531_199908-200910.nc', 'lon')/180*pi;
vwnd_all = ncread('vas_QuikSCAT_L2B_v20110531_199908-200910.nc', 'vas');

% extract the data from Jan. 2000 to Dec. 2008
uwnd_all = uwnd_all(:, :, 6:113);
vwnd_all = vwnd_all(:, :, 6:113);

% total number of months
T = 113-6+1;

n_lat = length(lat);
n_lon = length(lon);

[lat_m, lon_m] = meshgrid(lat, lon);

% the means of the u and v fields
m_uwnd = mean(uwnd_all, 3);
m_vwnd = mean(vwnd_all, 3);

% subtract the means
uwnd_resid = uwnd_all;
vwnd_resid = vwnd_all;
n_tot = n_lat*n_lon;
uwnd_resid_v = zeros(n_tot, T);
vwnd_resid_v = zeros(n_tot, T);
for i = 1:T
    uwnd_resid(:, :, i) = uwnd_resid(:, :, i)-m_uwnd;
    vwnd_resid(:, :, i) = vwnd_resid(:, :, i)-m_vwnd;
    uwnd_resid_v(:, i) = reshape(uwnd_resid(:, :, i), n_tot, 1);
    vwnd_resid_v(:, i) = reshape(vwnd_resid(:, :, i), n_tot, 1);
end

% reshape
m_uwnd_v = m_uwnd(:);
m_vwnd_v = m_vwnd(:);
lat_m_v = lat_m(:);
lon_m_v = lon_m(:);

% find the indices of NaN
idx_NaN = isnan(m_uwnd_v)+isnan(m_vwnd_v);

% extract the elements which are valid numbers
m_uwnd_v = m_uwnd_v(~idx_NaN);
m_vwnd_v = m_vwnd_v(~idx_NaN);
lat_m_v = lat_m_v(~idx_NaN);
lon_m_v = lon_m_v(~idx_NaN);
n_tot_val = length(m_uwnd_v);
uwnd_resid_val_v = zeros(n_tot_val, T);
vwnd_resid_val_v = zeros(n_tot_val, T);
for i = 1:T
    uwnd_resid_val_v(:, i) = uwnd_resid_v(~idx_NaN, i);
    vwnd_resid_val_v(:, i) = vwnd_resid_v(~idx_NaN, i);
end

% construct data matrix
X = [uwnd_resid_val_v' vwnd_resid_val_v'];
[U, S, V] = svd(X, 'econ');

% find the 95% threshold
lambda = diag(S).^2;
var_exp = cumsum(lambda)./sum(lambda);
thres = find(var_exp>=0.95, 1, 'first'); 

% compute the residual fields
V_L = V(:, 1:thres);
T_L = X*V_L;
X_L = T_L*V_L';
r = X-X_L;

% find the range of the region of Indian Ocean (IO)
% [1 2]*[-1 0.5]
range_lat_IO = find(lat>=-1, 1, 'first'):find(lat<=0.5, 1, 'last');
range_lon_IO = find(lon>=1, 1, 'first'):find(lon<=2, 1, 'last');

% plot the residual fields
h = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.075], [0.01 0.01], [0.075 0.01]);
i = 1;
u_field = r(i, 1:n_tot_val)';
v_field = r(i, n_tot_val+1:end)';
subplot(1, 2, 1)
u_field_m = NaN(n_tot, 1);
u_field_m(~idx_NaN) = u_field;
u_field_m = reshape(u_field_m, n_lon, n_lat);
u_field_m_IO = u_field_m(range_lon_IO, range_lat_IO);
lon_m_IO = lon_m(range_lon_IO, range_lat_IO);
lat_m_IO = lat_m(range_lon_IO, range_lat_IO);
pcolor(lon_m_IO/pi*180, lat_m_IO/pi*180, u_field_m_IO)
shading flat
axis equal
axis tight
colorbar
colormap(hot)
c_range = max(abs(min(u_field_m_IO(:))), max(u_field_m_IO(:)));
caxis([-c_range, c_range])
title('Jan. 2000 U Residual Field [m/s]')
xlabel('East longitude')
ylabel('Latitude')
set(gca, 'FontSize', 12)
subplot(1, 2, 2)
v_field_m = NaN(n_tot, 1);
v_field_m(~idx_NaN) = v_field;
v_field_m = reshape(v_field_m, n_lon, n_lat);
v_field_m_IO = v_field_m(range_lon_IO, range_lat_IO);
pcolor(lon_m_IO/pi*180, lat_m_IO/pi*180, v_field_m_IO)
shading flat
axis equal
axis tight
colorbar
colormap(hot)
c_range = max(abs(min(v_field_m_IO(:))), max(v_field_m_IO(:)));
caxis([-c_range, c_range])
title('Jan. 2000 V Residual Field [m/s]')
xlabel('East longitude')
ylabel('Latitude')
set(gca, 'FontSize', 12)
set(h, 'Position', [0, 0, 550, 350]);

% save the locations
u_field_m_IO_v = u_field_m_IO(:);
v_field_m_IO_v = v_field_m_IO(:);
idx_NaN_IO = isnan(u_field_m_IO_v)+isnan(v_field_m_IO_v);
lon_m_IO_v = lon_m_IO(~idx_NaN_IO);
lat_m_IO_v = lat_m_IO(~idx_NaN_IO);
u_field_m_IO_v = u_field_m_IO(~idx_NaN_IO);
v_field_m_IO_v = v_field_m_IO(~idx_NaN_IO);
theta = pi/2-lat_m_IO_v;
phi = lon_m_IO_v;

% convert coordinates
[x, y, z] = trans_coord(theta, phi);

n = length(x);
save('loc.mat', 'x', 'y', 'z', 'n', 'theta', 'phi')

% downsampling with half resolution
range_lat = range_lat_IO(1:2:end);
range_lon = range_lon_IO(1:2:end);
lon_m_sub = lon_m(range_lon, range_lat);
lat_m_sub = lat_m(range_lon, range_lat);

% record the residual fields
for i = 1:T
    u_field = r(i, 1:n_tot_val)';
    v_field = r(i, n_tot_val+1:end)';
    u_field_m = NaN(n_tot, 1);
    u_field_m(~idx_NaN) = u_field;
    u_field_m = reshape(u_field_m, n_lon, n_lat);
    v_field_m = NaN(n_tot, 1);
    v_field_m(~idx_NaN) = v_field;
    v_field_m = reshape(v_field_m, n_lon, n_lat);
    u_field_m_sub = u_field_m(range_lon, range_lat);
    v_field_m_sub = v_field_m(range_lon, range_lat);
    u_field_m_sub_v = u_field_m_sub(:);
    v_field_m_sub_v = v_field_m_sub(:);
  
    if i==1
        idx_NaN_sub = isnan(u_field_m_sub_v)+isnan(v_field_m_sub_v);
        lon_m_sub_v = lon_m_sub(:);
        lat_m_sub_v = lat_m_sub(:);
        lon_m_sub_v = lon_m_sub_v(~idx_NaN_sub);
        lat_m_sub_v = lat_m_sub_v(~idx_NaN_sub);
        n_sub_val = length(lon_m_sub_v);
        rec_u_field_m_sub_v = zeros(n_sub_val, T);
        rec_v_field_m_sub_v = zeros(n_sub_val, T);
    end
    
    u_field_m_sub_v = u_field_m_sub_v(~idx_NaN_sub);
    v_field_m_sub_v = v_field_m_sub_v(~idx_NaN_sub);
    rec_u_field_m_sub_v(:, i) = u_field_m_sub_v;
    rec_v_field_m_sub_v(:, i) = v_field_m_sub_v;
end

% convert to co-latitude and longitude
theta = pi/2-lat_m_sub_v;
phi = lon_m_sub_v;

% convert coordinates
[x, y, z] = trans_coord(theta, phi);

n = length(x);
samples = zeros(T, 2*n);

for i = 1:T
    samples(i, :) = reshape([rec_u_field_m_sub_v(:, i) rec_v_field_m_sub_v(:, i)]', 1, 2*n);
end

% save as .mat file
save('wind.mat', 'x', 'y', 'z', 'n', 'samples', 'theta', 'phi')
