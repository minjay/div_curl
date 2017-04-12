% combine the results of bootstrap for NMG
clear

load('boot_NMG1.mat')
rec_beta_hat_all = rec_beta_hat;
load('boot_NMG2.mat')
rec_beta_hat_all = rec_beta_hat_all+rec_beta_hat;
load('boot_NMG3.mat')
rec_beta_hat_all = rec_beta_hat_all+rec_beta_hat;
load('boot_NMG11.mat')
rec_beta_hat_all = rec_beta_hat_all+rec_beta_hat;
load('boot_NMG22.mat')
rec_beta_hat_all = rec_beta_hat_all+rec_beta_hat;

rec_beta_hat_all(:, 6) = 1./rec_beta_hat_all(:, 6);

boot_ste = std(rec_beta_hat_all, 0, 1);