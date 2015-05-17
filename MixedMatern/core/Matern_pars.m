function mat = Matern_pars(r, sigma1, sigma2, rho12, nu1, nu2, a)
%MATERN_PARS   Evaluates the parsimonious bivariate Matern cross-
%covariance function.
%
%   MAT = MATERN_PARS(R, SIGMA1, SIGMA2, RHO12, NU1, NU2, A);
%
% Inputs:
%   R - the Euclidean distance between two locations
%   SIGMA1, SIGMA2 - the standard deviation parameters
%   RHO12 - the co-located correlation
%   NU1, NU2 - the smoothness parameters
%   A - the spatial scale parameter
%
% Outputs:
%   MAT - the evaluated parsimonious bivariate Matern cross-covariance
%   function
%
% Author: Minjie Fan, 2015

mat = zeros(2);
mat(1, 1) = sigma1^2*Matern(r, nu1, a);
mat(2, 2) = sigma2^2*Matern(r, nu2, a);
mat(1, 2) = rho12*sigma1*sigma2*Matern(r, (nu1+nu2)/2, a);
mat(2, 1) = mat(1, 2);

end
