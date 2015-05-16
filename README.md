# div_curl
The codes accompanying the paper "Cross-covariance functions for divergence-free and curl-free tangent vector fields on the sphere"

## Dependencies
The codes are written in Matlab R2014a with Financial Toolbox installed. They also depend on the following packages written by other researchers.

* RBFSPHERE Matlab package http://math.boisestate.edu/~wright/montestigliano/index.html
* MEALPix 3.0 (The original link is broken)
* subtightplot http://www.mathworks.com/matlabcentral/fileexchange/39664-subtightplot
* Binned Scatter Plot http://www.mathworks.com/matlabcentral/fileexchange/19506-binned-scatter-plot

*To ensure compatibility, the latter three packages are already included in this repository.*

## Usage
1. Download the codes by cloning the repository `git clone https://github.com/minjay/div_curl`.
2. Install BFSPHERE Matlab package.

## Reproducing all numerical results in the paper
All scripts are in the folder `/MixedMatern/tests`, and all intermediate and final results are saved in the folder `/results`. Note that some codes take a **long** time to run even they are run parallelly on a server.

### Simulation part
1. To generate Figure 1: run `sim1.m` in the Matlab Command Window.
2. To generate Figure 2: first run `sim3_server.m` and `sim3_server2.m` on a server to obtain the MLEs of 500 realizations. The former uses a regular grid with (n_lat, n_lon)=(10, 20), while the latter uses a regular grid with (n_lat, n_lon)=(15, 30). After that, run `plot_sim.m` to generate Figure 2.
3. To generate Table 1: first run `sim4_bootstrap.m` to obtain the MLEs of 200 bootstrap samples. Then run `compute_sim_ratio.m` to compute the ratios between bootstrap standard errors and empirical standard errors.

### Real data part
1. The real data are downloaded through FTP: ftp://podaac-ftp.jpl.nasa.gov/allData/quikscat/L3/wind_1deg_1mo. The filenames are `uas_QuikSCAT_L2B_v20110531_199908-200910.nc` and `vas_QuikSCAT_L2B_v20110531_199908-200910.nc`. They represent u (zonal) and v (meridional) wind fields.
2. `real_init.m` preprocesses the real data, transforms the raw data to u and v residual fields, and generates Figure 3.
3. `real_sim.m` checks the multivariate normality of the u and v residual fields, and generates Figure 4.
4. `real_fit.m` fits the residual fields to the Mixed Matern model, and the results are summarized in the second column of Table 2, where the standard errors are computed by `real_bootstrap.m`. It takes 
5. 

