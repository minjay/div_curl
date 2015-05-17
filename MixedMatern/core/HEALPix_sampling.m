function [theta, phi, n] = HEALPix_sampling(Nside_power)
%HEALPIX_SAMPLING   Samples the observations on the HEALPix grid.
%
%   [THETA, PHI, N] = HEALPIX_SAMPLING(NSIDE_POWER);
%
% Inputs:
%   NSIDE_POWER - the number of grid points is NSIDE = 2^NSIDE_POWER
%
% Outputs:
%   THETA, PHI - the locations of the grid points in the spherical
%   coordinate
%   N - the number of grid points
%
% Author: Minjie Fan, 2015

Nside = 2^Nside_power;
tp = pix2ang(Nside, 'nest', false);
n = 12*Nside^2;
theta = zeros(n, 1);
phi = zeros(n, 1);
for i = 1:n
    theta(i) = tp{i}(1);
    phi(i) = tp{i}(2);
end

end

