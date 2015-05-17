% simulate div-free and curl-free random fields

clear

% HEALPix grid
[theta, phi, n] = HEALPix_sampling(3);

% convert coordinates
[x, y, z] = trans_coord(theta, phi);

% initial computation
[h_mat, r, P_cell, Q_cell, A_cell] = init_comp(x, y, z, n, theta, phi);

% generate samples
p = 2;
% set the same seed
rng('default')
samples = mvnrnd(zeros(p*n, 1), eye(p*n))';

% specify parameters
% [sigma1, sigma2, rho12, nu1, nu2, a, tau1, tau2]
beta_all = [1 0 0 2 3 3 0 0];
[u_curl1, v_curl1] = sim_mix(n, beta_all, samples, h_mat, r, P_cell, Q_cell, A_cell);

beta_all = [1 0 0 3 3 3 0 0];
[u_curl2, v_curl2] = sim_mix(n, beta_all, samples, h_mat, r, P_cell, Q_cell, A_cell);

beta_all = [0 1 0 3 2 3 0 0];
[u_div1, v_div1] = sim_mix(n, beta_all, samples, h_mat, r, P_cell, Q_cell, A_cell);

beta_all = [0 1 0 3 3 3 0 0];
[u_div2, v_div2] = sim_mix(n, beta_all, samples, h_mat, r, P_cell, Q_cell, A_cell);

% plot
h = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.01 0.01], [0.05 0.01]);
subplot(2, 2, 1)
plot_quivers(theta, phi, u_div1, v_div1)
title('Divergence-free with $$(\nu, a)=(2, 3)$$', 'interpreter', 'latex')
subplot(2, 2, 2)
plot_quivers(theta, phi, u_div2, v_div2)
title('Divergence-free with $$(\nu, a)=(3, 3)$$', 'interpreter', 'latex')
subplot(2, 2, 3)
plot_quivers(theta, phi, u_curl1, v_curl1)
title('Curl-free with $$(\nu, a)=(2, 3)$$', 'interpreter', 'latex')
subplot(2, 2, 4)
plot_quivers(theta, phi, u_curl2, v_curl2)
title('Curl-free with $$(\nu, a)=(3, 3)$$', 'interpreter', 'latex')
set(h, 'Position', [0, 0, 1000, 600]);
