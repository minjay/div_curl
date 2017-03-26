function value = fun_M(x, nu)
% compute function M_nu at x

if x>0
    value = x^nu*besselk(nu, x);
else
    value = gamma(nu)/(2^(1-nu));
end

end
