function mat = A_mat(theta, phi)

mat = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta);
    cos(theta)*cos(phi) cos(theta)*sin(phi) -sin(theta);
    -sin(phi) cos(phi) 0];

end