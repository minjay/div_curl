function c = get_c(h_mat, r, P_cell, Q_cell, A_cell, Matern_func, beta, coef, bessel, n_lat, n_lon, tau1, tau2)

c = cell(n_lat, n_lat, 3);

for i_lat = 1:n_lat
    for j_lat = 1:n_lat
        c_vec = zeros(n_lon, 3);
        for i_lon = 1:n_lon
            i_obs = (i_lat-1)*n_lon+i_lon;
            j_obs = (j_lat-1)*n_lon+1;
            if i_obs==j_obs
                mat = Matern_func(P_cell{j_obs}, Q_cell{j_obs}, P_cell{j_obs}, Q_cell{j_obs},...
                    A_cell{j_obs}, A_cell{j_obs}, h_mat{j_obs, j_obs}, r(j_obs, j_obs), beta, zeros(1, 3), zeros(1, 6));
                mat = (mat+mat')/2;
                mat = mat+diag([tau1^2, tau2^2]);
            else
                if i_obs<j_obs
                    h_mat_tmp = h_mat{i_obs, j_obs};
                    r_tmp = r(i_obs, j_obs);
                    idx = (j_obs-1)*(j_obs-2)/2+i_obs;
                else
                    h_mat_tmp = h_mat{j_obs, i_obs};
                    r_tmp = r(j_obs, i_obs);
                    idx = (i_obs-1)*(i_obs-2)/2+j_obs;
                end
                mat = Matern_func(P_cell{i_obs}, Q_cell{i_obs}, P_cell{j_obs}, Q_cell{j_obs},...
                    A_cell{i_obs}, A_cell{j_obs}, h_mat_tmp, r_tmp, beta, coef, bessel(idx, :));
            end
            c_vec(i_lon, 1) = mat(1, 1);
            c_vec(i_lon, 2) = mat(1, 2);
            c_vec(i_lon, 3) = mat(2, 2);
        end
        for l = 1:3
            c{i_lat, j_lat, l} = c_vec(:, l);
        end
    end
end

end

        