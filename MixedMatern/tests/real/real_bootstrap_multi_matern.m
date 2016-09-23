clear

load('wind.mat')

savefile = 'samples_multi_matern.mat';

T = 108;
p = 2;
% number of bootstrap samples
B = 200;

% initial computation
[h_mat, r, P_cell, Q_cell, A_cell] = init_comp(x, y, z, n, theta, phi);

param_BM = [0.21810846 0.203180249 0.058196367 1.23929461 1.1322247 0.184217146 0.156812964 -0.07977680];

beta_all(1:2) = param_BM(6:7);
beta_all(3) = param_BM(8);
beta_all(4:5) = param_BM(4:5);
beta_all(6) = param_BM(3);
beta_all(7:8) = param_BM(1:2);
sigma1 = sqrt(beta_all(1));
sigma2 = sqrt(beta_all(2));
rho12 = beta_all(3);
nu1 = beta_all(4);
nu2 = beta_all(5);
a = 1/beta_all(6);
tau1 = beta_all(7);
tau2 = beta_all(8);

% get cov mat
cov_mat = get_cov_Matern_pars(r, sigma1, sigma2, rho12, nu1, nu2, a)+...
    diag(kron(ones(1, n), [tau1^2, tau2^2]));

samples_all = mvnrnd(zeros(p*n, 1), cov_mat, B*T);

save(savefile, 'samples_all')
