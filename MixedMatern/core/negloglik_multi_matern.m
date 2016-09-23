function f_value = negloglik_multi_matern(beta_all, r, samples)
%NEGLOGLIK   A wrapper of ECMNOBJ to evaluate negative log-likelihood function. 
%
%   F_VALUE = NEGLOGLIK(BETA_ALL, R, SAMPLES);
%
% Inputs:
%   BETA_ALL - the current parameter vector
%   R - one of the outputs of the function INIT_COMP, which are computed 
%   beforehand to speed up the evaluation of the negative log-likelihood 
%   function
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
nu = (nu1+nu2)/2;
d = 3;
ub = gamma(nu1+d/2)^0.5*gamma(nu2+d/2)^0.5*gamma(nu)/gamma(nu1)^0.5/...
    gamma(nu2)^0.5/gamma(nu+d/2);
if abs(rho12)>ub
    beta_all(3) = sign(rho12)*ub;
end

disp(['Current estimate of beta is ', mat2str(round(beta_all*1e6)/1e6)])

sigma1 = beta_all(1);
sigma2 = beta_all(2);
a = beta_all(6);
tau1 = beta_all(7);
tau2 = beta_all(8);

% get cov mat
cov_mat = get_cov_Matern_pars(r, sigma1, sigma2, rho12, nu1, nu2, a)+...
    diag(kron(ones(1, n), [tau1^2, tau2^2]));
    
% negative log-likelihood function
f_value = ecmnobj(samples, zeros(p*n, 1), cov_mat);

end
