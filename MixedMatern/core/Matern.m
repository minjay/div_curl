function value = Matern(r, nu, a)
%MATERN   Evaluates the Matern correlation function.
%
%   VALUE = MATERN(R, NU, A);
%
% Inputs:
%   R - the Euclidean distance between two locations
%   NU - the smoothness parameter
%   A - the spatial scale parameter
%
% Outputs:
%   VALUE - the value of the Matern correlation function
%
% Author: Minjie Fan, 2015

if r>0
    value = 2^(1-nu)/gamma(nu)*(a*r)^nu*besselk(nu, a*r);
else
    value = 1;
end

end