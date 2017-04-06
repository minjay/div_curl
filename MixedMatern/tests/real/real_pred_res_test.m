% t-tests for prediction results

load('pred_err.mat')
res_TMM = [MSPE_u MSPE_v MAE_u MAE_v LogS_u LogS_v CRPS_u CRPS_v];
load('pred_err_BM.mat')
res_BM = [MSPE_u MSPE_v MAE_u MAE_v LogS_u LogS_v CRPS_u CRPS_v];
load('pred_err_NMG.mat')
res_NMG = [MSPE_u MSPE_v MAE_u MAE_v LogS_u LogS_v CRPS_u CRPS_v];

% compare TMM and NMG
n = size(res_TMM, 2);
p_TMM_NMG = zeros(n, 1);
for i = 1:n
    [h, p] = ttest2(res_TMM(:, i), res_NMG(:, i));
    p_TMM_NMG(i) = p;
end

% compare TMM and PARS-BM
p_TMM_BM = zeros(n, 1);
for i = 1:n
    [h, p] = ttest2(res_TMM(:, i), res_BM(:, i));
    p_TMM_BM(i) = p;
end
