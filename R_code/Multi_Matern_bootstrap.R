## This script is still under testing. It is problematic.

rm(list=ls())

library(RandomFields)
library(R.matlab)

RFoptions(modus_operandi="easygoing")
B <- 200

## read data
## change the directory if necessary
platform <- Sys.info()['sysname']
if (platform=='Darwin'){
  dat <- readMat('/Users/minjay/Documents/Documents/Research/div_curl/MixedMatern/tests/real/wind.mat')
  dat.boot <- readMat('/Users/minjay/Documents/Documents/Research/div_curl/results/samples_multi_matern.mat')
}else if (platform=='Linux'){
  dat <- readMat('/home/minjay/div_curl/MixedMatern/tests/real/wind.mat')
  dat.boot <- readMat('/home/minjay/div_curl/results/samples_multi_matern.mat')
}

samples <- dat$samples
samples.boot <- dat.boot$samples.all

n.rep <- dim(samples)[1]
n.obs <- dim(samples)[2]/2

## get distance matrix
## chordal distance
loc <- cbind(dat$x, dat$y, dat$z)
Dist.mat <- as.vector(dist(loc))

#################################
## model definition            ##
#################################
## bivariate pure nugget effect:
nug.sim <- RMmatrix(M = matrix(nc = 2, c(0.21810846, 0, 0, 0.203180249)), RMnugget())
## parsimonious bivariate Matern model
pars.model.sim <- nug.sim + RMbiwm(nudiag = c(1.23929461, 1.1322247), scale = 0.058196367, cdiag = c(0.184217146, 0.156812964),
                           rhored = -0.07977680)

## bivariate pure nugget effect:
nug <- RMmatrix(M = matrix(nc = 2, c(NA, 0, 0, NA)), RMnugget())
## parsimonious bivariate Matern model
pars.model <- nug + RMbiwm(nudiag = c(NA, NA), scale = NA, cdiag = c(NA, NA),
                           rhored = NA)

#################################
## parametric bootstrap        ##
#################################
for (i in 1:B){
  index <- ((i-1)*n.rep+1):(i*n.rep)
  samples <- samples.boot[index, ]
  U <- t(samples[, seq(1, n.obs*2, 2)])
  V <- t(samples[, seq(2, n.obs*2, 2)])
  UV <- matrix(0, n.obs, n.rep*2)
  UV[, seq(1, n.rep*2, 2)] <- U
  UV[, seq(2, n.rep*2, 2)] <- V
  # users.guess takes a model object
  pars <- RFfit(pars.model, distances = Dist.mat, dim = 3, data = UV, users.guess = pars.model.sim,
                optim.control = list(trace = TRUE))
}
