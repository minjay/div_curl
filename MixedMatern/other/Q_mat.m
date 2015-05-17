function mat = Q_mat(s)

mat = zeros(3);
mat(1, 2) = -s(3);
mat(1, 3) = s(2);
mat(2, 1) = s(3);
mat(2, 3) = -s(1);
mat(3, 1) = -s(2);
mat(3, 2) = s(1);

end