function [beta_hat, f_min] = Matern_fit_long(negloglik1, beta_init, lb, ub, mycon, paral)
%MATERN_FIT   Fits the Mixed Matern model by the MLE.
%
%   [BETA_HAT, F_MIN] = MATERN_FIT(NEGLOGLIK1, BETA_INIT, LB, UB, MYCON,
%   PARAL)
%
% Inputs:
%   NEGLOGLIK1 - the handle of the function evaluating the negative
%   log-likelihood
%   BETA_INIT - the initial guess of the parameter vector
%   LB, UB - the lower and upper bounds of the parameter vector
%   MYCON - additional nonlinear constraints of the parameter vector
%   PARAL - a boolean indicating whether gradients are estimated
%   parallelly in each iteration
%
% Outputs:
%   BETA_HAT - the optimal parameter vector minimizing the negative
%   log-likelihood
%   F_MIN - the value of the negative log-likelihood at BETA_HAT
%
% Author: Minjie Fan, 2015

disp(['The initial guess of beta is ', mat2str(round(beta_init*1e6)/1e6)])

% use interior-point algorithm for large-scale problems
% nu can not be too large

if ~paral
    options = optimoptions(@fmincon, 'Algorithm', 'interior-point',...
        'Hessian', 'lbfgs', 'Display', 'iter', 'TolX', 1e-6,...
        'MaxFunEvals', 20000);
else
    options = optimoptions(@fmincon, 'Algorithm', 'interior-point',...
        'Hessian', 'lbfgs', 'Display', 'iter', 'TolX', 1e-6,...
        'UseParallel', 'always', 'MaxFunEvals', 20000);
end

[beta_hat, f_min] = fmincon(negloglik1, beta_init, [], [], [], [], lb, ub,...
    mycon, options);

disp(['The MLE of beta is ', mat2str(round(beta_hat*1e6)/1e6)])

end