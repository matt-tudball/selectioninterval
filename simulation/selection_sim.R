rm(list=ls())
require(ggplot2); require(parallel); require(devtools); require(nloptr);
require(numDeriv); require(pbapply); library(survey)
setwd("C:/Users/ow18301/OneDrive - University of Bristol/MyFiles-Migrated/Documents/IEASBOS_FILES_")

# Load functions
source("CODE/selectioninterval/R/selection_bound.R")
source("CODE/selectioninterval/R/misc_funs.R")

# Seed
set.seed(2021)

# Population size
N <- 1e7

# Population data
X <- as.matrix(rbeta(N,2,5))
Y <- as.vector(runif(N,0,1))
W <- as.matrix(cbind(X,Y))

cons <- list(
  list('RESP', 0.07),
  list('COVMEAN', 1, 0.5)#,
  #list('DIREC', 1, '+')
)

# Population confidence interval
results_pop <- selection_bound(y=Y, x=X, w=W, L0l=0.1, L0u=0.2, L1=3)

# Simulation
sample <- c(seq(10,50,5), seq(100,1000,100), seq(1500,10000,500))
main <- rep(NA, length(sample))

for (j in 1:length(sample)) {
  if (get_os() == "windows") {
    cl <- makeCluster(4, outfile = "")
    clusterExport(cl, c("selection_bound","X","Y","W","N","sample","j","auglag","grad",
                        "jacobian","results_pop","nloptr","svydesign","svyglm"),
                  envir = environment())
  } else {
    cl <- nnodes
  }

  out <- t(pbsapply(X=1:1e4, cl = cl, simplify=T, FUN=function(x) {
    # Sample data
    n <- sample[j]
    id <- sample(N, n, replace=T, prob=NULL)
    Ys <- Y[id]; Xs <- X[id,]; Ws <- W[id,]

    # Compute sample confidence interval
    results_sample <- selection_bound(y=Ys, x=Xs, w=Ws, L0l=0.1, L0u=0.2, L1=3)

    # Coverage
    coverage <- (results_sample$ci[1] <= results_pop$interval[1]) &
      (results_sample$ci[2] >= results_pop$interval[2])

    return(coverage)
  }))

  if (get_os() == "windows") {
    stopCluster(cl)
  }

  main[j] <- mean(out)
  print(paste("Coverage is ", mean(out)," for sample size ", sample[j], sep=""))
}

# Make the plot
data <- data.frame(x=sample, y=main)
plot <- ggplot(data=data, aes(x=x, y=y)) +
  geom_line() +
  geom_point() +
  xlab("Sample size") +
  ylab("Coverage frequency") +
  scale_y_continuous(breaks=seq(0,1,0.2)) +
  geom_hline(yintercept=0.95, color="grey42")
#ggsave(filename="FAMMR_FILES/FIGURES/power_curves.pdf",plot=plot)

