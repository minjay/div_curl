% compute the stderr ratios

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

var_names = {'$$\sigma_{1}$$', '$$\sigma_{2}$$', '$$\rho_{12}$$', '$$\nu_{1}$$', '$$\nu_{2}$$', '$$a$$', '$$\tau_{1}$$', '$$\tau_{2}$$'};
boxplot(ratio, 'Labels', var_names)
ylabel('bootstrap stderr/empirical stderr')
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'FontSize', 12);

ratio_median = median(ratio, 1);