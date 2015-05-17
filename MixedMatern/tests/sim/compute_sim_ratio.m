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

ratio_median = median(ratio, 1);