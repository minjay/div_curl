function [h_mat, r, P_cell, Q_cell, A_cell] = init_comp(x, y, z, n, theta, phi)

h_mat = cell(n);
r = zeros(n);
% column first
for j = 1:n
    for i = 1:j
        s = [x(i); y(i); z(i)];
        t = [x(j); y(j); z(j)];
        h_mat{i, j} = (s-t)*(s-t)';
        r(i ,j) = norm(s-t);
    end
end

P_cell = cell(n, 1);
Q_cell = cell(n, 1);
for i = 1:n
    s = [x(i); y(i); z(i)];
    P_cell{i} = P_mat(s);
    Q_cell{i} = Q_mat(s);
end

A_cell = cell(n, 1);
for i = 1:n
    A_cell{i} = [0 0 1;0 -1 0]*A_mat(theta(i), phi(i));
end

end
