rm(list=ls())
library(devtools); library(foreign)
setwd("C:/Users/ow18301/OneDrive - University of Bristol/MyFiles-Migrated/Documents/IEASBOS_FILES_")

# Load function
load_all(path='CODE/selectioninterval')

# Load data
dat <- read.dta("DATA/data_for_example.dta")

# Format data
x <- dat[,colnames(dat)[startsWith(colnames(dat),'x')]]
z <- dat[,c('z_rosla', colnames(x)[colnames(x) != 'x_edu'])]
w <- dat[,colnames(dat)[startsWith(colnames(dat),'w')]]
y <- dat[,'y_inc']

# Constraints
mycons <- list(
  list('RESP', 0.055),
  list('COVMEAN', 2, 0.495),
  list('DIREC', 1, '+')
)

out1 <- selection_bound(y=y, x=x, z=z, w=w, L0l=0.02, L0u=0.2, L1=2, cons=NULL, theta=NULL, alpha=0.05)

out2 <- selection_bound(y=y, x=x, z=z, w=w, L0l=0.02, L0u=0.2, L1=2, cons=mycons, theta=NULL, alpha=0.05)
