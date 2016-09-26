load('MLE1.mat')
rec_beta_hat1 = rec_beta_hat;
load('MLE2.mat')
rec_beta_hat2 = rec_beta_hat;
load('MLE7.mat')
rec_beta_hat3 = rec_beta_hat;

h = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0.075], [0.1 0.05], [0.05 0.02]);

var_names = {'\sigma_{1}', '\sigma_{2}', '\rho_{12}', '\nu_{1}', '\nu_{2}', '1/a', '\tau_{1}', '\tau_{2}'};
lab_cell = {'10*20', '15*30', '20*40'};
true_values = [1 1 0.5 3 4 1/2 0.1 0.1];
for i = 1:8
    subplot(2,4,i)
    if i==6
        myboxplot([1./rec_beta_hat1(:,i) 1./rec_beta_hat2(:,i) 1./rec_beta_hat3(:,i)], lab_cell)
    else
        myboxplot([rec_beta_hat1(:,i) rec_beta_hat2(:,i) rec_beta_hat3(:,i)], lab_cell)
    end
    title(['$$', var_names{i}, '$$'], 'interpreter', 'latex')
    line([0.5 4.5], [true_values(i) true_values(i)], 'Color', 'k', 'LineStyle', '--')
    hline = findobj(gca, 'tag', 'Median');
    set(hline, 'LineWidth', 1.25)
    set(gca, 'FontSize', 12)
end

set(h, 'Position', [0, 0, 800, 500]);
