clear

load('wind.mat')

savefile = 'samples_NMG.mat';

T = 108;
p = 2;
% number of bootstrap samples
B = 200;

% initial computation
[r, h0_cell] = init_comp_NMG(n, theta, phi, x, y, z);

beta_all = [0.016283 -0.004629 0.002649 0.011662 2.867271 14.523555 0.270156 0.343079 0.720569 0.830202 0.212515 0.194964];
a1 = beta_all(1);
a2 = beta_all(2);
b1 = beta_all(3);
b2 = beta_all(4);
nu = beta_all(5);
a = beta_all(6);
sigma1 = beta_all(7);
sigma2 = beta_all(8);
w1 = beta_all(9);
w2 = beta_all(10);
tau1 = beta_all(11);
tau2 = beta_all(12);
cov_mat = get_cov_NMG(r, a1, a2, b1, b2, nu, a, h0_cell)+...
    get_cov_Matern_pars(r, sigma1, sigma2, 0, w1, w2, a)+...
    diag(kron(ones(1, n), [tau1^2, tau2^2]));

samples_all = mvnrnd(zeros(p*n, 1), cov_mat, B*T);

samples_all_cell = cell(B, 1);
for rep = 1:B
    samples_all_cell{rep} = samples_all((rep-1)*T+1:rep*T, :);
end

save(savefile, 'samples_all_cell');
