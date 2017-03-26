function value = fun_M(x, nu)
% compute function M_nu at x

if x>0
    value = x^nu*besselk(nu, x);
elseif nu~=0
    value = gamma(nu)/(2^(1-nu));
else
    value = 0;
end

end
