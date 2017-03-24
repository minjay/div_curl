function lh = boxplot_curve(r_upper, cov_u_upper, nb, col, y_range)

[X_MED,Y_MED,Y_LOW,Y_HIGH] = binned_plot(r_upper, cov_u_upper, nb, 'y_range', y_range);

hold on
for i = 1:nb
    if i==1
        lh = line([X_MED(i), X_MED(i)], [Y_LOW(i), Y_HIGH(i)], 'Color', col, 'LineWidth', 1);
    else
        line([X_MED(i), X_MED(i)], [Y_LOW(i), Y_HIGH(i)], 'Color', col, 'LineWidth', 1)
    end
end

end