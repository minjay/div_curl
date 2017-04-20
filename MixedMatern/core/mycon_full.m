function [c, ceq] = mycon_full(beta_all)

rho12 = beta_all(3);
nu1 = beta_all(4);
nu2 = beta_all(5);
nu12 = beta_all(6);
d = 3;
ub = gamma(nu1+d/2)^0.5*gamma(nu2+d/2)^0.5*gamma(nu12)/gamma(nu1)^0.5/...
    gamma(nu2)^0.5/gamma(nu12+d/2);
c = [rho12-ub -rho12-ub (nu1+nu2)/2-nu12];
ceq = [];

end