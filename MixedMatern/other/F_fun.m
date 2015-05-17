function value = F_fun(r, nu, a, coef, bessel)

if r>0
    value = -coef*bessel;
else
    value = -a^2/(nu-1)/2;
end

end