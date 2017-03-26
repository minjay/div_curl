% run on server
parpool(8)
addpath(genpath('/home/minjay/div_curl'))

load('wind.mat')

% initial computation
r = zeros(n);
h10 = zeros(n);
h20 = zeros(n);
h30 = zeros(n);
h120 = zeros(n);
h130 = zeros(n);
h230 = zeros(n);
h330 = zeros(n);
% column first
for j = 1:n
    for i = 1:j
        s = [x(i); y(i); z(i)];
        t = [x(j); y(j); z(j)];
        r(i, j) = norm(s-t);
        [h10_tmp, h20_tmp, h30_tmp, h120_tmp, h130_tmp, h230_tmp, h330_tmp] = h_deriv(1, theta(i), theta(j), phi(i)-phi(j));
        h10(i, j) = h10_tmp;
        h20(i, j) = h20_tmp;
        h30(i, j) = h30_tmp;
        h120(i, j) = h120_tmp;
        h130(i, j) = h130_tmp;
        h230(i, j) = h230_tmp;
        h330(i, j) = h330_tmp;
    end
end

% negative log-likelihood function
negloglik1 = @(beta_all) negloglik_NMG_Matern(beta_all, r, samples,...
    h10, h20, h30, h120, h130, h230, h330);

lb = [0 -10 -10 -10 1 10 0 0 1 1 0 0];
ub = [10 10 10 10 5 10 10 10 5 5 10 10];
beta_init = rand_search_NMG(negloglik1, 10, lb, ub, true, true);

% to avoid identifiability problem, set a1>0
lb = [0 -Inf -Inf -Inf 1 0 0 0 1 1 0 0];
ub = [Inf Inf Inf Inf 5 Inf Inf Inf 5 5 Inf Inf];

% fit the model
[beta_hat, f_min] = Matern_fit(negloglik1, beta_init, lb, ub, [], true);

delete(gcp)
