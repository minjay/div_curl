function [phi_diff, cov_data_sub, cov_sub_TMM, cov_sub_BM, cov_sub_NMG,...
    lat1, lat2] = extract_cov(i1, i2, cov_data, cov_TMM, cov_BM, cov_NMG,...
    lat, lon)
% function to extract a subset of covariance matrix

lat_unique = unique(lat);

lat1 = lat_unique(i1);
lat2 = lat_unique(i2);
index1 = find(lat==lat1);
index2 = find(lat==lat2);

cov_data_sub = [];
cov_sub_TMM = [];
cov_sub_BM = [];
cov_sub_NMG = [];
phi_diff = [];

for i = 1:length(index1)
    for j = 1:length(index2)
        ii = index1(i);
        jj = index2(j);
        % extract the pairs on the specified latitudes
        cov_data_sub = [cov_data_sub cov_data(ii, jj)];
        cov_sub_TMM = [cov_sub_TMM cov_TMM(ii, jj)];
        cov_sub_BM = [cov_sub_BM cov_BM(ii, jj)];
        cov_sub_NMG = [cov_sub_NMG cov_NMG(ii, jj)];
        phi_diff = [phi_diff lon(ii)-lon(jj)];
    end
end

end