function dist = chordal_dist(theta1, theta2, phi)
% compute chordal dist

theta = theta1-theta2;
dist = 2*sqrt((sin(theta/2))^2+sin(theta1)*sin(theta2)*(sin(phi/2))^2);

end