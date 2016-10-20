function cov_mat = get_cov_full(h_mat, r, P_cell, Q_cell, A_cell, Matern_func, beta, coef, bessel)
%GET_COV   Gets the cross-covariance matrix of the Mixed Matern model
%
%   COV_MAT = GET_COV(H_MAT, R, P_CELL, Q_CELL, A_CELL, MATERN_FUNC, BETA,
%   COEF, BESSEL);
%
% Inputs:
%   H_MAT, R, P_CELL, Q_CELL, A_CELL - the outputs of the function
%   INIT_COMP, which are computed beforehand to speed up the computation of
%   the cross-covariance matrix
%   MATERN_FUNC - the handle of the function to evaluate the Mixed Matern
%   cross-covariance function
%   BETA - the current parameter vector (excluding tau_1 and tau_2)
%   COEF, BESSEL - the outputs of the function GET_COEF_BESSEL, which are
%   computed beforehand to speed up the computation of the cross-covariance
%   matrix
%
% Outputs:
%   COV_MAT - the obtained cross-covariance matrix
%
% Author: Minjie Fan, 2015

n = size(r, 1);
p = 2;
cov_mat = zeros(p*n);
idx = 0;

% column first
for j = 1:n
    for i = 1:j-1
        idx = idx+1;
        mat = Matern_func(P_cell{i}, Q_cell{i}, P_cell{j}, Q_cell{j}, A_cell{i}, A_cell{j}, h_mat{i, j}, r(i, j), beta, coef, bessel(idx, :));
        cov_mat(p*(i-1)+1:p*i, p*(j-1)+1:p*j) = mat;
        mat = Matern_func(P_cell{j}, Q_cell{j}, P_cell{i}, Q_cell{i}, A_cell{j}, A_cell{i}, h_mat{j, i}, r(j, i), beta, coef, bessel(idx, :));
        cov_mat(p*(j-1)+1:p*j, p*(i-1)+1:p*i) = mat;
    end      
    % when i=j
    mat = Matern_func(P_cell{j}, Q_cell{j}, P_cell{j}, Q_cell{j}, A_cell{j}, A_cell{j}, h_mat{j, j}, r(j, j), beta, zeros(1, 3), zeros(1, 6));
    % make it sym
    mat = (mat+mat')/2;
    cov_mat(p*(j-1)+1:p*j, p*(j-1)+1:p*j) = mat;      
end

end
