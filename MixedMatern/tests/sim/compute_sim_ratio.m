% compute the stderr ratios
clear

stderr_boot = zeros(10, 8);
load('MLE2.mat')
stderr_emp = std(rec_beta_hat);

for i = 1:10
    load(['boot_sim', num2str(i), '.mat'])
    stderr_boot(i, :) = std(rec_beta_hat);
end

ratio = zeros(10, 8);
for i = 1:10
    ratio(i, :) = stderr_boot(i, :)./stderr_emp;
end

subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0.075], [0.05 0.03], [0.05 0.02]);
subplot(1, 1, 1)
var_names = {'$$\sigma_{1}$$', '$$\sigma_{2}$$', '$$\rho_{12}$$', '$$\nu_{1}$$', '$$\nu_{2}$$', '$$a$$', '$$\tau_{1}$$', '$$\tau_{2}$$'};
boxplot(ratio, 'Labels', var_names)
ylabel('Ratio')
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'FontSize', 12);

ratio_median = median(ratio, 1);