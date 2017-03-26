function f_value = negloglik_NMG(beta_all, r, samples, h10,...
    h20, h30, h120, h130, h230, h330)

n = size(r, 1);
p = 2;

disp(['Current estimate of beta is ', mat2str(round(beta_all*1e6)/1e6)])

a1 = beta_all(1);
a2 = beta_all(2);
b1 = beta_all(3);
b2 = beta_all(4);
nu = beta_all(5);
a = beta_all(6);
tau1 = beta_all(7);
tau2 = beta_all(8);

% get cov mat
cov_mat = get_cov_NMG(r, a1, a2, b1, b2, nu, a, h10, h20,...
    h30, h120, h130, h230, h330)+...
    diag(kron(ones(1, n), [tau1^2, tau2^2]));
    
% negative log-likelihood function
f_value = ecmnobj(samples, zeros(p*n, 1), cov_mat);

end
