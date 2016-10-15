% simulate the data for examining the prediction performance

clear

% HEALPix grid
[theta, phi, n] = HEALPix_sampling(3);

% convert coordinates
[x, y, z] = trans_coord(theta, phi);

% initial computation
[h_mat, r, P_cell, Q_cell, A_cell] = init_comp(x, y, z, n, theta, phi);

p = 2;
% set the same seed
rng('default')

% specify parameters
% [sigma1, sigma2, rho12, nu1, nu2, a, tau1, tau2]
beta_all = [1 1 0.5 3 4 2 0.1 0.1];

samples = mvnrnd(zeros(p*n, 1), eye(p*n))';
[u, v] = sim_mix(n, beta_all, samples, h_mat, r, P_cell, Q_cell, A_cell);
samples = reshape([u v]', 1, p*n);

% CV locations
B = 500;
width = pi/6;
n_pred = n/2;
n_est = n-n_pred;
rec_pred_loc = zeros(B, n_pred);
rec_est_loc = zeros(B, n_est);
rec_idx_est = zeros(B, p*n_est);
rec_idx_pred = zeros(B, p*n_pred);
for i = 1:B
    lb = rand*2*pi;
    ub = lb+width;
    if ub>2*pi
        ub_adj = ub-2*pi;
        pop = find(phi>ub_adj & phi<lb);
    else
        pop = find(phi<lb | phi>ub);
    end
    rec_est_loc(i, :) = sort(randsample(pop, n_est));
    rec_pred_loc(i, :) = setdiff(1:n, rec_est_loc(i, :));

    rec_idx_est(i, :) = sort([rec_est_loc(i, :)*p-1 rec_est_loc(i, :)*p]);
    rec_idx_pred(i, :) = setdiff(1:p*n, rec_idx_est(i, :));
end

filename = 'sim_data_mix.mat';
save(filename, 'samples', 'theta', 'phi', 'n', 'x', 'y', 'z',...
    'rec_pred_loc', 'rec_est_loc', 'rec_idx_est', 'rec_idx_pred')
