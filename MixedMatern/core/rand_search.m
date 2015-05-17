function beta_init = rand_search(negloglik1, n_cand, lb, ub, lhs, paral)
%RAND_SEARCH   Searches for the best initial value for the MLE.
%
%   BETA_INIT = RAND_SEARCH(NEGLOGLIK1, N_CAND, LB, UB, LHS, PARAL);
%
% Inputs:
%   NEGLOGLIK1 - the handle of the function used to evaluate the negative
%   log-likelihood
%   N_CAND - the number of candidates
%   LB, UB - the lower and upper bounds of the parameter vector
%   LHS - whether use Latin hypercube sampling or uniform random sampling
%   PARAL - whether the function is run parallelly or not
%
% Outputs:
%   BETA_INIT - the best initial parameter vector
%
% Author: Minjie Fan, 2015

n_para = length(lb);
beta_cand = zeros(n_para, n_cand);
if ~lhs
    for i = 1:n_para
        beta_cand(i, :) = rand(1, n_cand)*(ub(i)-lb(i))+lb(i);
    end
else
    X = lhsdesign(n_cand, n_para);
    for i = 1:n_para
        beta_cand(i, :) = X(:, i)'*(ub(i)-lb(i))+lb(i);
    end
end

d = 3;

% make sure that pho12 satisfies the constraint
for i = 1:n_cand
    nu1 = beta_cand(4, i);
    nu2 = beta_cand(5, i);
    nu = (nu1+nu2)/2;
    ub = gamma(nu1+d/2)^0.5*gamma(nu2+d/2)^0.5*gamma(nu)/gamma(nu1)^0.5/...
        gamma(nu2)^0.5/gamma(nu+d/2);
    if abs(beta_cand(3, i))>ub
        beta_cand(3, i) = sign(beta_cand(3, i))*ub;
    end
end

tmp_value = zeros(1, n_cand);
if paral
    parfor i = 1:n_cand
        tmp_value(i) = negloglik1(beta_cand(:, i));
    end
else
    for i = 1:n_cand
        tmp_value(i) = negloglik1(beta_cand(:, i));
    end
end

[~, index] = min(tmp_value);
beta_init = beta_cand(:, index(1))';

end