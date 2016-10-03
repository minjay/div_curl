function [] = myboxplot(rec_beta_hat, lab_cell, loc)

if nargin<3
    loc = 0.15;
end

boxplot(rec_beta_hat, 'labels', lab_cell);
[~, c] = size(rec_beta_hat);
med = median(rec_beta_hat);
hold on
for i = 1:c
    textString = sprintf('%.2f', med(i));
    text(i-loc, med(i), textString);
end

end