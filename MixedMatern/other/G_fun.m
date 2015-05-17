function value = G_fun(r, coef, bessel)

if r>0
    value = coef*bessel;
else
    value = 0;
end

end