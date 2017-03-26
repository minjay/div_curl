function value = fun_M(x, nu)
% compute function M_nu at x

value = x^nu*besselk(nu, x);

end
