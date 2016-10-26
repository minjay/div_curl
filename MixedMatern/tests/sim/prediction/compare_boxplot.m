clear

load('sim_pred_rep_TMM.mat')

MSPE_u_TMM = MSPE_u;
MSPE_v_TMM = MSPE_v;
MAE_u_TMM = MAE_u;
MAE_v_TMM = MAE_v;
LogS_u_TMM = LogS_u;
LogS_v_TMM = LogS_v;
CRPS_u_TMM = CRPS_u;
CRPS_v_TMM = CRPS_v;

load('pred_err_BM_sim.mat')

subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0.1], [0.05 0.05], [0.05 0.05]);

models = {'TMM', 'PARS_BM'};

subplot(2, 4, 1)
boxplot([MSPE_u_TMM MSPE_u], 'labels', models)
title('MSPE, u')
subplot(2, 4, 2)
boxplot([MSPE_v_TMM MSPE_v], 'labels', models)
title('MSPE, v')
subplot(2, 4, 3)
boxplot([MAE_u_TMM MAE_u], 'labels', models)
title('MAE, u')
subplot(2, 4, 4)
boxplot([MAE_v_TMM MAE_v], 'labels', models)
title('MAE, v')
subplot(2, 4, 5)
boxplot([LogS_u_TMM LogS_u], 'labels', models)
title('LogS, u')
subplot(2, 4, 6)
boxplot([LogS_v_TMM LogS_v], 'labels', models)
title('LogS, v')
subplot(2, 4, 7)
boxplot([CRPS_u_TMM CRPS_u], 'labels', models)
title('CRPS, u')
subplot(2, 4, 8)
boxplot([CRPS_v_TMM CRPS_v], 'labels', models)
title('CRPS, v')
