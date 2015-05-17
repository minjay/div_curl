function [x, y, z] = trans_coord(theta, phi)
%TRANS_COORD   Transforms the spherical coordinate to the Euclidean
%coordinate.
%
%   [X, Y, Z] = TRANS_COORD(THETA, PHI);
%
% Inputs:
%   THETA, Phi - the spherical coordinate
%
% Outputs:
%   X, Y, Z - the Euclidean coordinate
%
% Author: Minjie Fan, 2015

x = sin(theta).*cos(phi);
y = sin(theta).*sin(phi);
z = cos(theta);

end