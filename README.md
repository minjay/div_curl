# div_curl
The codes accompanying the paper "Modeling tangential vector fields on a sphere". **Note that this description is out of date and needs to be updated.**

The aim of this repository is to allow the reproduction of all numerical results reported in this paper.

## Dependencies
Except the analyses by the multivariate Matern model (based on the R package RandomFields), all the codes are written in Matlab R2014a with Financial Toolbox installed. They also depend on the following packages written by other researchers.

* RBFSPHERE Matlab package http://math.boisestate.edu/~wright/montestigliano/index.html
* MEALPix 3.0 (The original link is broken)
* subtightplot http://www.mathworks.com/matlabcentral/fileexchange/39664-subtightplot
* Binned Scatter Plot http://www.mathworks.com/matlabcentral/fileexchange/19506-binned-scatter-plot

*To ensure compatibility, the latter three packages are already included in the repository.*

## Usage
1. Download the codes by cloning the repository `git clone https://github.com/minjay/div_curl`.
2. Install BFSPHERE Matlab package by following the instructions on http://math.boisestate.edu/~wright/montestigliano/index.html.

## Reproducing all numerical results in the paper
All main scripts are in the folders `/MixedMatern/tests` (Matlab) and `/R_code` (R), and all intermediate and final results are saved in the folder `/results`. Note that some codes take a **long** time to run even they are run parallelly on a server.

### Simulation part
1. To generate Figure 1: run `sim1.m` in the Matlab Command Window.
2. To generate Figure 2: first run `sim3_server.m` and `sim3_server2.m` on a server to obtain the MLEs of 500 realizations. The former uses a regular grid with (n_lat, n_lon)=(10, 20), while the latter uses a regular grid with (n_lat, n_lon)=(15, 30). After that, run `plot_sim.m` to generate Figure 2.
3. To generate Table 1: first, run the command `matlab -r "sim4_bootstrap(num)"` on a server to obtain the MLEs of 200 bootstrap samples, where `num` is 1 to 10 since this procedure is done for 10 realizations. Then run `compute_sim_ratio.m` to compute the ratios between bootstrap standard errors and empirical standard errors and the median ratios for these 10 realizations.

### Real data part
1. The real data are downloaded from the FTP: ftp://podaac-ftp.jpl.nasa.gov/allData/quikscat/L3/wind_1deg_1mo. The filenames are `uas_QuikSCAT_L2B_v20110531_199908-200910.nc` and `vas_QuikSCAT_L2B_v20110531_199908-200910.nc`. They represent u (zonal) and v (meridional) wind fields, respectively.
2. `real_init.m` preprocesses the real data, transforms the raw data to u and v residual fields, and generates Figure 3.
3. `real_sim.m` checks the multivariate normality of the u and v residual fields, and generates Figure 4.
4. `real_fit.m` fits the residual fields to the Mixed Matern model, and the results are summarized in the second column of Table 2, where the standard errors are computed by `real_bootstrap.m`. Both of them are run on a server. Whether the physical characteristics of the wind fields are captured by the Mixed Matern model is checked by `real_sim.m`. It generates Figure 5 showing the empirical and fitted cross-covariance functions.
5. The residual fields are also fitted to the parsimonious bivariate Matern model by `Multi_Matern.R`. It uses the R package RandomFields and runs on a server. The results are summarized in the first column of Table 2.
6. Figure 6 is plotted by `real_kriging.m`.
7. The cokriging cross-validation results shown in Table 3 are generated by `real_kriging20.m` (for the Mixed Matern model) and `Multi_Matern_kriging.R`, `real_kriging20_BM.m` (for the parsimonious bivariate Matern model). The first two scripts are run on a server.

Some results would vary slightly because different random numbers could be generated.

