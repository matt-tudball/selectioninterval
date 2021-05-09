rm(list=ls())
require(ggplot2); require(mvtnorm); require(devtools); require(nloptr); require(numDeriv)
setwd("C:/Users/ow18301/OneDrive - University of Bristol/MyFiles-Migrated/Documents/IEASBOS_FILES_")

# Load functions
source("CODE/selectioninterval/R/opt_funs.R")
source("CODE/selectioninterval/R/selection_bound.R")

# Population size
n <- 1000

# Observed data
X <- rnorm(n,0,1)
Y <- 1 + 0.5*X + rnorm(n,0,sqrt(1-0.5^2))
W <- cbind(X,Y)

cons <- list(
  list('RESP', 0.3),
  list('COVMEAN', X, -0.5)
)

# Confidence interval
results <- selection_bound(y=Y, x=X, w=W, L0=3, L1=2, cons=cons)

print(results$value)

theta1 <- results$par
inv_wgt <- 1+exp(-W%*%theta1)
mean(inv_wgt)
mean(inv_wgt*X)/mean(inv_wgt)

results <- auglag(x0=c(0,0,0), fn=fn, gr=gr, lower=c(log(1/L0), rep(log(1/L1),p)),
                  upper=c(log(L0), rep(log(L1),p)), hin=NULL, hinjac=NULL,
                  localsolver=c("SLSQP"), nl.info=F, control=nlopts)

print(results$value)
