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
1. To generate Figure 1: run `sim1.m` in the Matlab Command Window.
2. To generate Figure 2: first run `sim3_server.m` and `sim3_server2.m` on a server to obtain the MLEs of 500 realizations. The former uses a regular grid with (n_lat, n_lon)=(10, 20), while the latter with (n_lat, n_lon)=(15, 30). Note that it takes a long time to finish running. Then 


