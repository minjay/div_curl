function [theta, phi, n] = regular_sampling(n_lat, n_lon)
%REGULAR_SAMPLING   Samples the observations on a regular grid of the
%   sphere. The latitudes range from -50 degree to 50 degree.
%
%   [THETA, PHI, N] = REGULAR_SAMPLING(N_LON, N_LAT);
%
% Inputs:
%   N_LON - the number of longitudes
%   N_LAT - the number of latitudes
%
% Outputs:
%   THETA, PHI - the locations of the grid points in the spherical
%   coordinate
%   N - the number of grid points
%
% Author: Minjie Fan, 2015

st = 50/180*pi;
en = 130/180*pi;
d_lat = (en-st)/(n_lat+1);
d_lon = 2*pi/n_lon;
theta = linspace(st+d_lat, en-d_lat, n_lat);
phi = linspace(0, 2*pi-d_lon, n_lon);
[theta, phi] = meshgrid(theta, phi);
theta = theta(:);
phi = phi(:);
n = length(theta);

end