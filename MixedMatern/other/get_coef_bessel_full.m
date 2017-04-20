function [coef, bessel] = get_coef_bessel(beta, r)

nu1 = beta(4);
nu2 = beta(5);
nu = (nu1+nu2)/2;
a = beta(6);
nu_vec = [nu1 nu2 nu];
coef = 2.^(1-nu_vec)./gamma(nu_vec)*a^2;
r_vec = r(r~=0);
bessel = zeros(length(r_vec), 6);
bessel(:, 1) = (a*r_vec).^(nu1-1).*besselk(nu1-1, a*r_vec);
bessel(:, 2) = (a*r_vec).^(nu2-1).*besselk(nu2-1, a*r_vec);
bessel(:, 3) = (a*r_vec).^(nu-1).*besselk(nu-1, a*r_vec);
bessel(:, 4) = a^2*(a*r_vec).^(nu1-2).*besselk(nu1-2, a*r_vec);
bessel(:, 5) = a^2*(a*r_vec).^(nu2-2).*besselk(nu2-2, a*r_vec);
bessel(:, 6) = a^2*(a*r_vec).^(nu-2).*besselk(nu-2, a*r_vec);

end