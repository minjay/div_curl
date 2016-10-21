function mat_div = Matern_div(P_mat_s, Q_mat_s, P_mat_t, Q_mat_t, A_mat_s, A_mat_t, h_mat, r, beta, coef, bessel)
%MATERN_DIV   Evaluates Cov(div-free(s), y(t)).
%
%   MAT_DIV = MATERN_DIV(P_MAT_S, Q_MAT_S, P_MAT_T, Q_MAT_T, A_MAT_S, A_MAT_T,
%   H_MAT, R, BETA, COEF, BESSEL);
%
% Inputs:
%   P_MAT_S, Q_MAT_S, P_MAT_T, Q_MAT_T, A_MAT_S, A_MAT_T - the values of
%   the matrices P, Q, A at s and t
%   H_MAT - the matrix of (s-t)*(s-t)'
%   R - the Euclidean distance between s and t
%   BETA - the current parameter vector (excluding tau_1 and tau_2)
%   COEF, BESSEL - two quantities computed beforehand by the function 
%   GET_COEF_BESSEL to speed up the evaluation of the Mixed Matern 
%   cross-covariance function
%
% Outputs:
%   MAT_DIV
%
% Author: Minjie Fan, 2015

sigma1 = beta(1);
sigma2 = beta(2);
rho12 = beta(3);
nu1 = beta(4);
nu2 = beta(5);
a = beta(6);
nu = (nu1+nu2)/2;

% only include Qs*Pt' and Qs*Qt' terms
mat_div = -sigma2^2*Q_mat_s*K_fun(h_mat, r, nu2, a, coef(2), bessel([2 5]))*Q_mat_t'-...
    rho12*sigma1*sigma2*Q_mat_s*K_fun(h_mat, r, nu, a, coef(3), bessel([3 6]))*P_mat_t';
mat_div = A_mat_s*mat_div*A_mat_t';

end