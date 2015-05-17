function value = K_fun(h_mat, r, nu, a, coef, bessel)

value = F_fun(r, nu, a, coef, bessel(1))*eye(3)+G_fun(r, coef, bessel(2))*h_mat;

end