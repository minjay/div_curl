function cov_mat = get_cov_Matern_pars(r, sigma1, sigma2, rho12, nu1, nu2, a)
%GET_COV_MATERN_PARS   Gets the cross-covariance matrix of the parsimonious
%bivariate Matern model
%
%   COV_MAT = GET_COV_MATERN_PARS(R, SIGMA1, SIGMA2, RHO12, NU1, NU2, A);
%
% Inputs:
%   R - one of the outputs of the function INIT_COMP, which records the 
%   Euclidean distances between all the pairs of sampling locations
%   SIGMA1, SIGMA2 - the standard deviation parameters
%   RHO12 - the co-located correlation
%   NU1, NU2 - the smoothness parameter
%   A - the spatial scale parameter
%
% Outputs:
%   COV_MAT - the obtained cross-covariance matrix
%
% Author: Minjie Fan, 2015

n = size(r, 1);
p = 2;
cov_mat = zeros(p*n);

for j = 1:n
    for i = 1:j
        mat = Matern_pars(r(i, j), sigma1, sigma2, rho12, nu1, nu2, a);
        cov_mat(p*(i-1)+1:p*i, p*(j-1)+1:p*j) = mat;
        cov_mat(p*(j-1)+1:p*j, p*(i-1)+1:p*i) = mat';
    end
end

end
