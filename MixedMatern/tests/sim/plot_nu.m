% verify why nu2 is estimated much worse than nu1

clear

load('MLE1.mat')
rec_beta_hat1 = rec_beta_hat;
% increase sigma2
load('MLE10.mat')
rec_beta_hat10 = rec_beta_hat;
% decrease nu2
load('MLE11.mat')
rec_beta_hat11 = rec_beta_hat;

subplot(1, 2, 1)
boxplot([rec_beta_hat1(:, 4) rec_beta_hat11(:, 4) rec_beta_hat10(:, 4) ], 'labels',...
    {'$$\sigma_2=1, \nu_2=4$$', '$$\sigma_2=\sqrt{2/3}, \nu_2=3$$', '$$\sigma_2=2, \nu_2=4$$'});
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
title('\nu_1')
subplot(1, 2, 2)
boxplot([rec_beta_hat1(:, 5) rec_beta_hat11(:, 5) rec_beta_hat10(:, 5)], 'labels',...
    {'$$\sigma_2=1, \nu_2=4$$', '$$\sigma_2=\sqrt{2/3}, \nu_2=3$$', '$$\sigma_2=2, \nu_2=4$$'});
title('\nu_2')
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);