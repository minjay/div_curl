function cov_mat = get_cov_NMG(r, a1, a2, b1, b2, nu, a, h0_cell)

n = size(r, 1);
p = 2;
cov_mat = zeros(p*n);

for j = 1:n
    for i = 1:j
        mat = NMG_fast(r(i, j), [a1 a2], [b1 b2], nu, a, h0_cell{i, j});
        cov_mat(p*(i-1)+1:p*i, p*(j-1)+1:p*j) = mat;
        cov_mat(p*(j-1)+1:p*j, p*(i-1)+1:p*i) = mat';
    end
end

end
