function [r, h0_cell] = init_comp_NMG(n, theta, phi, x, y, z)

% initial computation
r = zeros(n);
% an n-by-n cell
h0_cell = cell(n);
% column first
for j = 1:n
    for i = 1:j
        s = [x(i); y(i); z(i)];
        t = [x(j); y(j); z(j)];
        % chordal distance
        r(i, j) = norm(s-t);
        % obtain h's when a=1
        [h1, h2, h3, h12, h13, h23, h33] = h_deriv(1, theta(i), theta(j), phi(i)-phi(j));
        h_vec = [h1, h2, h3, h12, h13, h23, h33];
        h0_cell{i, j} = h_vec;
    end
end

end