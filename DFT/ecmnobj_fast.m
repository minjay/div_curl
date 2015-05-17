function f_value = ecmnobj_fast(samples, c, n_lat, n_lon)
%ECMNOBJ_FAST   A fast version of ECMNOBJ to evaluate negative
%log-likelihood function by the DFT.
%
%   F_VALUE = ECMNOBJ_FAST(SAMPLES, C, N_LAT, N_LON);
%
% Inputs:
%   SAMPLES - the row vector of stacked (u, v) observations with length
%   p*n, where p=2
%   C - the cell recording the elements of the cross-covariance matrix,
%   which can be obtained by the function GET_C
%   N_LAT - the number of latitudes
%   N_LON - the number of longitudes
%
% Outputs:
%   F_VALUE - the value of the negative log-likelihood function
%
% Author: Minjie Fan, 2015

p = 2;
u = reshape(samples(1:p:end), n_lon, n_lat);
v = reshape(samples(2:p:end), n_lon, n_lat);

% transform u and v
trans_u = trans_y(u);
trans_v = trans_y(v);
trans_u = trans_u(:);
trans_v = trans_v(:);
new_u = zeros(n_lat, n_lon);
new_v = zeros(n_lat, n_lon);

for i_lat = 1:n_lat
    for k = 1:n_lon
        new_u(i_lat, k) = trans_u(n_lon*(i_lat-1)+k);
        new_v(i_lat, k) = trans_v(n_lon*(i_lat-1)+k);
    end
end

blk_diag_mat = cell(n_lon, 3);
for l = 1:3
    for k = 1:n_lon
        blk_diag_mat{k, l} = zeros(n_lat);
    end
end

for i_lat = 1:n_lat
    for j_lat = 1:n_lat
        for l = 1:3
            diag_elem = fft(c{i_lat, j_lat, l});
            % re-arrange
            for k = 1:n_lon
                blk_diag_mat{k, l}(i_lat, j_lat) = diag_elem(k);
            end
        end
    end
end

D1 = blk_diag_mat(:, 1);
D12 = blk_diag_mat(:, 2);
D2 = blk_diag_mat(:, 3);

inv_D1 = inv_blk_diag(D1);
inv_D2 = inv_blk_diag(D2);
log_det_D2 = log_det_blk_diag(D2);
H_D12 = H_blk_diag(D12);

mix_mat1 = sub_blk_diag(D1, multi_blk_diag(multi_blk_diag(D12, inv_D2), H_D12));
mix_mat2 = sub_blk_diag(D2, multi_blk_diag(multi_blk_diag(H_D12, inv_D1), D12));

log_det_mix_mat1 = log_det_blk_diag(mix_mat1);
log_det_Sigma = log_det_mix_mat1+log_det_D2;

A = inv_blk_diag(mix_mat1);
D = inv_blk_diag(mix_mat2);
B = gmultiply(-1, multi_blk_diag(multi_blk_diag(inv_D1, D12), D));
C = H_blk_diag(B);

tmp = 0;
for k = 1:n_lon
    tmp = tmp+new_u(:, k)'*A{k}*new_u(:, k)+new_v(:, k)'*C{k}*new_u(:, k)+...
        new_u(:, k)'*B{k}*new_v(:, k)+new_v(:, k)'*D{k}*new_v(:, k);
end

f_value = real((log_det_Sigma+tmp+n_lon*n_lat*p*log(2.0*pi))/2.0);

end

