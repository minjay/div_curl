function [] = myboxplot(rec_beta_hat, lab_cell)

boxplot(rec_beta_hat, 'labels', lab_cell);
[~, c] = size(rec_beta_hat);
med = median(rec_beta_hat);
hold on
for i = 1:c
    textString = sprintf('%.2f', med(i));
    text(i-0.15, med(i), textString);
end

end