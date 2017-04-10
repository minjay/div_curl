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

load('sim_pred_rep_NMG.mat')

MSPE_u_NMG = MSPE_u;
MSPE_v_NMG = MSPE_v;
MAE_u_NMG = MAE_u;
MAE_v_NMG = MAE_v;
LogS_u_NMG = LogS_u;
LogS_v_NMG = LogS_v;
CRPS_u_NMG = CRPS_u;
CRPS_v_NMG = CRPS_v;

load('pred_err_BM_sim.mat')

h = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.125 0.075], [0.05 0.05], [0.03 0.02]);

models = {'TMM', 'PARS_BM', 'NBG'};

subplot(2, 4, 1)
boxplot([MSPE_u_TMM MSPE_u MSPE_u_NMG], 'labels', models)
title('MSPE, u')
subplot(2, 4, 2)
boxplot([MSPE_v_TMM MSPE_v MSPE_v_NMG], 'labels', models)
title('MSPE, v')
subplot(2, 4, 3)
boxplot([MAE_u_TMM MAE_u MAE_u_NMG], 'labels', models)
title('MAE, u')
subplot(2, 4, 4)
boxplot([MAE_v_TMM MAE_v MAE_v_NMG], 'labels', models)
title('MAE, v')
subplot(2, 4, 5)
boxplot([LogS_u_TMM LogS_u LogS_u_NMG], 'labels', models)
title('LogS, u')
subplot(2, 4, 6)
boxplot([LogS_v_TMM LogS_v LogS_v_NMG], 'labels', models)
title('LogS, v')
subplot(2, 4, 7)
boxplot([CRPS_u_TMM CRPS_u CRPS_u_NMG], 'labels', models)
title('CRPS, u')
subplot(2, 4, 8)
boxplot([CRPS_v_TMM CRPS_v CRPS_v_NMG], 'labels', models)
title('CRPS, v')

set(h, 'Position', [0, 0, 700, 350]);
