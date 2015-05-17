function [u, v] = sim_mix(n, beta_all, samples, h_mat, r, P_cell, Q_cell, A_cell)

beta = beta_all(1:6);
tau1 = beta_all(7);
tau2 = beta_all(8);

% get coef and bessel
[coef, bessel] = get_coef_bessel(beta, r);

% get cov mat
cov_mat = get_cov(h_mat, r, P_cell, Q_cell, A_cell, @Matern_mix,...
    beta, coef, bessel)+diag(kron(ones(1, n), [tau1^2, tau2^2]));

L = chol(cov_mat, 'lower');

new_samples = L*samples;

p = 2;
u = new_samples(1:p:end);
v = new_samples(2:p:end);

end

