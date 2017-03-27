% sanity check for MLE fitting
rng('default')

[theta, phi, n] = HEALPix_sampling(3);

% convert coordinates
[x, y, z] = trans_coord(theta, phi);

[r, h0_cell] = init_comp_NMG(n, theta, phi, x, y, z);

p = 2;

% specify parameters
% [a1 a2 b1 b2 nu a tau1 tau2]
a1 = 1;
a2 = -1;
b1 = 1;
b2 = -1;
nu = 4;
a = 2;
tau1 = 0.1;
tau2 = 0.1;

beta_all = [a1 a2 b1 b2 nu a tau1 tau2];

cov_mat = get_cov_NMG(r, a1, a2, b1, b2, nu, a, h0_cell)+...
    diag(kron(ones(1, n), [tau1^2, tau2^2]));

samples = mvnrnd(zeros(1, p*n), cov_mat);

% negative log-likelihood function
negloglik1 = @(beta_all) negloglik_NMG(beta_all, r, samples, h0_cell);

negloglik1(beta_all)

beta_init = [0 0 0 0 3 3 0.2 0.2];
lb = [-2 -2 -2 -2 1 0 0 0];
ub = [2 2 2 2 5 10 10 10];
% fit the model
[beta_hat, f_min] = Matern_fit(negloglik1, beta_init, lb, ub, [], false);
