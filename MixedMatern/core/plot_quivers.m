function plot_quivers(theta, phi, u, v)
%PLOT_QUIVERS   Plots simulated samples on the sphere (after Hammer
%   projection).
%
%   PLOT_QUIVERS(THETA, PHI, U, V);
%
% Inputs:
%   THETA, PHI - the sampling locations
%   U, V - the u and v components
%
% Author: Minjie Fan, 2015

% convert theta and phi
theta = pi/2-theta;
phi(phi>pi)=phi(phi>pi)-2*pi;

[HX, HY] = sph2hammer(phi, theta);
quiver(HX, HY, u, v, 'b');

hold on

% draw contour
th = linspace(-pi/2,pi/2,101);
lam = -pi+0*th;
[xh,yh] = sph2hammer(lam,th);
plot(xh,yh,'k');
lam = pi+0*th;
[xh,yh] = sph2hammer(lam,th);
plot(xh,yh,'k');

axis equal
axis tight
axis off

set(gca, 'FontSize', 12)

end