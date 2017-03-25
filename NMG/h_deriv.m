function [h1, h2, h3, h12, h13, h23, h33] = h_deriv(a, theta1, theta2, phi)
% compute h1, h2, h3, h12, h13, h23, h33

theta = theta1-theta2;
two_a_sq = 2*a^2;
sin_phi_over_2_sq = (sin(phi/2))^2;
h1 = two_a_sq*(sin(theta)+2*cos(theta1)*sin(theta2)*sin_phi_over_2_sq);
h2 = two_a_sq*(-sin(theta)+2*sin(theta1)*cos(theta2)*sin_phi_over_2_sq);
h3 = two_a_sq*sin(theta1)*sin(theta2)*sin(phi);
h12 = two_a_sq*(-cos(theta)+2*cos(theta1)*cos(theta2)*sin_phi_over_2_sq);
h13 = two_a_sq*cos(theta1)*sin(theta2)*sin(phi);
h23 = two_a_sq*sin(theta1)*cos(theta2)*sin(phi);
h33 = two_a_sq*sin(theta1)*sin(theta2)*cos(phi);

end