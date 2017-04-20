function [c, ceq] = mycon(beta_all)

rho12 = beta_all(3);
nu1 = beta_all(4);
nu2 = beta_all(5);
nu = (nu1+nu2)/2;
d = 3;
ub = gamma(nu1+d/2)^0.5*gamma(nu2+d/2)^0.5*gamma(nu)/gamma(nu1)^0.5/...
    gamma(nu2)^0.5/gamma(nu+d/2);
c = [rho12-ub -rho12-ub];
ceq = [];

end