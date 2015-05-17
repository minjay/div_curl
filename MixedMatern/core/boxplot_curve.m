function boxplot_curve(r_upper, cov_u_upper, nb, col)

[X_MED,Y_MED,Y_LOW,Y_HIGH] = binned_plot(r_upper, cov_u_upper, nb);

hold on
for i = 1:nb
    line([X_MED(i), X_MED(i)], [Y_LOW(i), Y_HIGH(i)], 'Color', col, 'LineWidth', 1)
end

end