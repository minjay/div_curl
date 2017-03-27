function mat = NMG_fast(r, a_vec, b_vec, nu, a, h0_vec)

%   R - the Euclidean distance between two locations (theta1, phi1) and
%   (theta2, phi2)

mat = zeros(2);
h_sqrt = a*r;
M_nu_minus_2 = fun_M(h_sqrt, nu-2);
M_nu_minus_1 = fun_M(h_sqrt, nu-1);

% h_vec = [h1 h2 h3 h12 h13 h23 h33];
% multiply a^2 to update since for h0_vec, a=1
h_vec = a^2*h0_vec;
h1 = h_vec(1);
h2 = h_vec(2);
h3 = h_vec(3);
h12 = h_vec(4);
h13 = h_vec(5);
h23 = h_vec(6);
h33 = h_vec(7);

for i = 1:2
    for j = 1:2
        ai = a_vec(i);
        aj = a_vec(j);
        bi = b_vec(i);
        bj = b_vec(j);
        % obtain Gamma coefficients
        [Gamma1, Gamma2] = coef_Gamma(ai, aj, bi, bj, h1, h2, h3, h12, h13, h23, h33);
        mat(i, j) = Gamma1*M_nu_minus_2+Gamma2*M_nu_minus_1;
    end
end

end
