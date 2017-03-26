% sanity checking for NMG

b_vec = [1 0.1];
% standardize by the Earth radius
a = 1/(2000/6371);
nu = 5;
theta1 = pi/2;
phi1 = 0;
phi2 = pi/6;
phi = phi1-phi2;

theta2_vec = 0:1e-3:pi;
n = length(theta2_vec);
record_corr1 = zeros(n, 4);
record_corr2 = zeros(n, 4);
a2_all = [0 0.1 1 10];
for t = 1:4
    a_vec = [1 a2_all(t)];
    for i = 1:n
        theta2 = theta2_vec(i);
        r = chordal_dist(theta1, theta2, phi);
        mat = NMG(theta1, phi1, theta2, phi2, r, a_vec, b_vec, nu, a);
        % cov(Z1(theta1, phi1), Z1(theta1, phi1))
        mat2 = NMG(theta1, phi1, theta1, phi1, 0, [a_vec(1) a_vec(1)], [b_vec(1) b_vec(1)], nu, a);
        % cov(Z2(theta2, phi2), Z2(theta2, phi2))
        mat3 = NMG(theta2, phi2, theta2, phi2, 0, [a_vec(2) a_vec(2)], [b_vec(2) b_vec(2)], nu, a);
        % corr
        record_corr1(i, t) = mat(1, 2)/sqrt(mat2(1, 1)*mat3(1, 1));
        mat = NMG(theta2, phi2, theta1, phi1, r, a_vec, b_vec, nu, a);
        % cov(Z1(theta2, phi2), Z1(theta2, phi2))
        mat2 = NMG(theta2, phi2, theta2, phi2, 0, [a_vec(1) a_vec(1)], [b_vec(1) b_vec(1)], nu, a);
        % cov(Z2(theta1, phi1), Z2(theta1, phi1))
        mat3 = NMG(theta1, phi1, theta1, phi1, 0, [a_vec(2) a_vec(2)], [b_vec(2) b_vec(2)], nu, a);
        % corr
        record_corr2(i, t) = mat(1, 2)/sqrt(mat2(1, 1)*mat3(1, 1));
    end
end

plot(theta2_vec, record_corr1-record_corr2)
% adjust axis
axis([0 pi -0.5 0.5])

