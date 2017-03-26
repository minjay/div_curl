function cov_mat = get_cov_NMG(r, a1, a2, b1, b2, nu, a, h10, h20, h30,...
    h120, h130, h230, h330)

n = size(r, 1);
p = 2;
cov_mat = zeros(p*n);

for j = 1:n
    for i = 1:j
        mat = NMG_fast(r(i, j), [a1 a2], [b1 b2], nu, a, h10(i, j),...
            h20(i, j), h30(i, j), h120(i, j), h130(i, j), h230(i, j),...
            h330(i, j));
        cov_mat(p*(i-1)+1:p*i, p*(j-1)+1:p*j) = mat;
        cov_mat(p*(j-1)+1:p*j, p*(i-1)+1:p*i) = mat';
    end
end

end
