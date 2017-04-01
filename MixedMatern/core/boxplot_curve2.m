function lh = boxplot_curve2(r_upper, cov_u_upper, nb, col, y_range, width)

[X_MED,Y_MED,Y_LOW,Y_HIGH] = binned_plot(r_upper, cov_u_upper, nb, 'y_range', y_range);

hold on
for i = 1:nb
    if i==1
        lh = line([X_MED(i)-width, X_MED(i)+width], [Y_LOW(i), Y_LOW(i)], 'Color', col, 'LineWidth', 1);
        line([X_MED(i)-width, X_MED(i)+width], [Y_HIGH(i), Y_HIGH(i)], 'Color', col, 'LineWidth', 1);
    else
        line([X_MED(i)-width, X_MED(i)+width], [Y_LOW(i), Y_LOW(i)], 'Color', col, 'LineWidth', 1);
        line([X_MED(i)-width, X_MED(i)+width], [Y_HIGH(i), Y_HIGH(i)], 'Color', col, 'LineWidth', 1);
    end
end

end