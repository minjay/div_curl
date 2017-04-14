function [ph1, ph2, ph3] = plot_three(phi_diff_sorted, res_TMM, res_BM, res_NMG)
% helper function for real_nonsta_axial_sym.m

ph1 = plot(phi_diff_sorted, res_TMM, 'b', 'LineWidth', 1.5);
ph2 = plot(phi_diff_sorted, res_BM, 'k-.', 'LineWidth', 1.5);
ph3 = plot(phi_diff_sorted, res_NMG, 'r--', 'LineWidth', 1.5);

end