function value = CRPS( x, mu, sigma_sq )
%CRPS Evaluates CRPS values

sigma = sqrt(sigma_sq);
x0 = (x-mu)./sigma;
cdf_x0 = normcdf(x0, 0, 1);
pdf_x0 = normpdf(x0, 0, 1);
value = sigma.*(x0.*(2*cdf_x0-1)+2*pdf_x0-1/sqrt(pi));

end

