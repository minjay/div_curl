% 15*30, tau=0.1
load('MLE2.mat')
rec_beta_hat1 = rec_beta_hat;
% tau=0.2
load('MLE3.mat')
rec_beta_hat2 = rec_beta_hat;
% tau=0.3
load('MLE4.mat')
rec_beta_hat3 = rec_beta_hat;

h = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0.075], [0.1 0.05], [0.05 0.02]);

var_names = {'\sigma_{1}', '\sigma_{2}', '\rho_{12}', '\nu_{1}', '\nu_{2}', '1/a', '\tau_{1}', '\tau_{2}'};
lab_cell = {'low', 'medium', 'high'};
true_values = [1 1 0.5 3 4 1/2];
for i = 1:8
    subplot(2,4,i)
    if i==6
        myboxplot([1./rec_beta_hat1(:,i) 1./rec_beta_hat2(:,i) 1./rec_beta_hat3(:,i)], lab_cell)
    else
        myboxplot([rec_beta_hat1(:,i) rec_beta_hat2(:,i) rec_beta_hat3(:,i)], lab_cell)
    end
    title(['$$', var_names{i}, '$$'], 'interpreter', 'latex')
    if i<=6
        line([0.5 4.5], [true_values(i) true_values(i)], 'Color', 'k', 'LineStyle', '--')
    else
        line([0.5 1.5], [0.1 0.1], 'Color', 'k', 'LineStyle', '--')
        line([1.5 2.5], [0.2 0.2], 'Color', 'k', 'LineStyle', '--')
        line([2.5 3.5], [0.3 0.3], 'Color', 'k', 'LineStyle', '--')
    end
    hline = findobj(gca, 'tag', 'Median');
    set(hline, 'LineWidth', 1.25)
    set(gca, 'FontSize', 12)
end

set(h, 'Position', [0, 0, 800, 500]);
