load('MLE_covariate.mat')
rec_beta_hat3 = rec_beta_hat;
load('MLE_covariate2.mat')
rec_beta_hat2 = rec_beta_hat;
load('MLE_covariate3.mat')
rec_beta_hat1 = rec_beta_hat;

subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0.075], [0.05 0.05], [0.03 0.02]);

var_names = {'\sigma_{1}', '\sigma_{2}', '\rho_{12}', '\nu_{1}', '\nu_{2}',...
    '1/a', '\tau_{1}', '\tau_{2}', '\beta_{10}', '\beta_{11}', '\beta_{12}',...
    '\beta_{20}', '\beta_{21}', '\beta_{22}'};

lab_cell = {'10*20', '15*30', '20*40'};
true_values = [1 1 0.5 3 4 1/2 0.1 0.1  5 5 -1 10 -15 3];
loc = [0.3 * ones(1, 10) 0.4 0.3 0.4 0.3];

h = figure;
for i = 1:14
    subplot(3, 5, i)
    if i==6
        myboxplot([1./rec_beta_hat1(:,i) 1./rec_beta_hat2(:,i) 1./rec_beta_hat3(:,i)], lab_cell, loc(i))
    else
        myboxplot([rec_beta_hat1(:, i) rec_beta_hat2(:, i) rec_beta_hat3(:, i)], lab_cell, loc(i))
    end
    title(['$$', var_names{i}, '$$'], 'interpreter', 'latex')
    line([0.5 4.5], [true_values(i) true_values(i)], 'Color', 'k', 'LineStyle', '--')
    hline = findobj(gca, 'tag', 'Median');
    set(hline, 'LineWidth', 1.25)
    set(gca, 'FontSize', 12)
    axis square
end

set(h, 'Position', [0, 0, 800, 500]);  