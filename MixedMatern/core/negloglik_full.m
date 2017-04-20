function f_value = negloglik_full(beta_all, h_mat, r, P_cell, Q_cell, A_cell, samples)
%NEGLOGLIK   A wrapper of ECMNOBJ to evaluate negative log-likelihood function. 
%
%   F_VALUE = NEGLOGLIK(BETA_ALL, H_MAT, R, P_CELL, Q_CELL, A_CELL,
%   SAMPLES);
%
% Inputs:
%   BETA_ALL - the current parameter vector
%   H_MAT, R, P_CELL, Q_CELL, A_CELL - the outputs of the function
%   INIT_COMP, which are computed beforehand to speed up the evaluation of
%   the negative log-likelihood function
%   SAMPLES - the T-by-n*p matrix of observations, where T is the number of
%   replicates, n is the number of sampling locations and p=2
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
nu12 = beta_all(6);
%a = beta_all(7);
d = 3;
ub = gamma(nu1+d/2)^0.5*gamma(nu2+d/2)^0.5*gamma(nu12)/gamma(nu1)^0.5/...
    gamma(nu2)^0.5/gamma(nu12+d/2);
if abs(rho12)>ub
    beta_all(3) = sign(rho12)*ub;
end

beta = beta_all(1:7);

disp(['Current estimate of beta is ', mat2str(round(beta_all*1e6)/1e6)])
[coef, bessel] = get_coef_bessel(beta, r);
tau1 = beta_all(8);
tau2 = beta_all(9);
cov_mat = get_cov(h_mat, r, P_cell, Q_cell, A_cell, @Matern_mix_full,...
    beta, coef, bessel)+diag(kron(ones(1, n), [tau1^2, tau2^2]));
    
% negative log-likelihood function
f_value = ecmnobj(samples, zeros(p*n, 1), cov_mat);

end
