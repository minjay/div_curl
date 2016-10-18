rm(list=ls())

library(RandomFields)
library(R.matlab)
require(parallel)

B <- 20

## two warnings
## seem to be bugs

## version of "RandomFields": 3.0.62

RFoptions(modus_operandi="sloppy")

## read data
## change the directory if necessary
platform <- Sys.info()['sysname']
if (platform=='Darwin'){
  dat <- readMat('/Users/minjay/Documents/MATLAB/Needlets/Vector/div_curl/MixedMatern/tests/real/wind.mat')
  pred.loc <- readMat('/Users/minjay/Documents/MATLAB/Needlets/Vector/div_curl/MixedMatern/tests/real/pred_loc.mat')
}else if (platform=='Linux'){
  dat <- readMat('/home/minjay/div_curl/MixedMatern/tests/real/wind.mat')
  pred.loc <- readMat('/home/minjay/div_curl/MixedMatern/tests/real/pred_loc.mat')
}

## re-arrange the data
samples <- dat$samples
n.rep <- dim(samples)[1]
n.obs <- dim(samples)[2]/2
U <- t(samples[, seq(1, n.obs*2, 2)])
V <- t(samples[, seq(2, n.obs*2, 2)])
UV <- matrix(0, n.obs, n.rep*2)
UV[, seq(1, n.rep*2, 2)] <- U
UV[, seq(2, n.rep*2, 2)] <- V

loc <- cbind(dat$x, dat$y, dat$z) 

rec.pred.loc <- pred.loc$rec.pred.loc
rec.est.loc <- pred.loc$rec.est.loc

#################################
## model definition            ##
#################################
## bivariate pure nugget effect:
nug <- RMmatrix(M = matrix(nc = 2, c(NA, 0, 0, NA)), RMnugget())
## parsimonious bivariate Matern model
pars.model <- nug + RMbiwm(nudiag = c(NA, NA), scale = NA, cdiag = c(NA, NA),
                           rhored = NA)
## whole bivariate Matern model
whole.model <- nug + RMbiwm(nudiag = c(NA, NA), nured = NA, s = rep(NA, 3),
                            cdiag = c(NA, NA), rhored = NA)

#################################
## model fitting and testing   ##
#################################
## parsimonious model
fit.model <- function(rep){
  pred.loc <- rec.pred.loc[rep, ]
  est.loc <- rec.est.loc[rep, ]
  
  loc.sub <- loc[est.loc, ]
  Dist.mat <- as.vector(dist(loc.sub))
  UV.sub <- UV[est.loc, ]
  
  pars <- RFfit(pars.model, distances = Dist.mat, dim = 3, data = UV.sub,
                optim.control = list(trace = TRUE))
  param <- print(pars)$param
  return(param[1, ])
}

param.BM <- mclapply(1:B, fit.model, mc.cores = 10)
param.BM <- t(matrix(unlist(param.BM), 8, B))
writeMat('param_kriging.mat', param_BM = param.BM)
