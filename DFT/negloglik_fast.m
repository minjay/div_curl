function f_value = negloglik_fast(beta_all, h_mat, r, P_cell, Q_cell, A_cell, samples, n_lat, n_lon)
%NEGLOGLIK_FAST   A wrapper of ECMNOBJ_FAST.
%
%   F_VALUE = NEGLOGLIK_FAST(BETA_ALL, H_MAT, R, P_CELL, Q_CELL, A_CELL, SAMPLES, N_LAT, N_LON)
%
% Inputs: 
%   BETA_ALL - the current parameter vector
%   H_MAT, R, P_CELL, Q_CELL, A_CELL - the outputs of the function
%   INIT_COMP, which are computed beforehand to speed up the evaluation of
%   negative log-likelihood function
%   SAMPLES - the row vector of stacked (u, v) observations with length
%   p*n, where p=2
%   n_lat - the number of latitudes
%   n_lon - the number of longitudes
%
% Outputs:
%   F_VALUE - the value of the negative log-likelihood function
%
% Author: Minjie Fan, 2015

n = size(r, 1);
p = 2;

rho12 = beta_all(3);
nu1 = beta_all(4);
nu2 = beta_all(5);
nu = (nu1+nu2)/2;
d = 3;
ub = gamma(nu1+d/2)^0.5*gamma(nu2+d/2)^0.5*gamma(nu)/gamma(nu1)^0.5/...
    gamma(nu2)^0.5/gamma(nu+d/2);
if abs(rho12)>ub
    beta_all(3) = sign(rho12)*ub;
end

beta = beta_all(1:6);

disp(['Current estimate of beta is ', mat2str(round(beta_all*1e6)/1e6)])
[coef, bessel] = get_coef_bessel(beta, r);
tau1 = beta_all(7);
tau2 = beta_all(8);

% get c
c = get_c(h_mat, r, P_cell, Q_cell, A_cell, @Matern_mix,...
    beta, coef, bessel, n_lat, n_lon, tau1, tau2);
    
% negative log-likelihood function
f_value = ecmnobj_fast(samples, c, n_lat, n_lon);

end