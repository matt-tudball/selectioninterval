rm(list=ls())
require(ggplot2); require(parallel); require(devtools); require(nloptr);
require(numDeriv); require(pbapply)
setwd("C:/Users/ow18301/OneDrive - University of Bristol/MyFiles-Migrated/Documents/IEASBOS_FILES_")

# Load functions
source("CODE/selectioninterval/R/selection_bound.R")
source("CODE/selectioninterval/R/misc_funs.R")

# Population size
N <- 1e6

# Population data
X <- as.matrix(rnorm(N,0,1))
Y <- as.vector(1 + 0.5*X + rnorm(N,0,sqrt(1-0.5^2)))
W <- as.matrix(cbind(X,Y))

cons <- list(
  list('RESP', 0.07),
  list('COVMEAN', X, 0.5)
)

# Population confidence interval
results_pop <- selection_bound(y=Y, x=X, w=W, L0l=0.05, L0u=0.15, L1=3)

# Simulation
sample <- 1000 #c(seq(10,100,5), seq(110,140,10), seq(150,1000,50), seq(1500,5000,500))
main <- rep(NA,length(sample))

for (j in 1:length(sample)) {
  if (get_os() == "windows") {
    cl <- makeCluster(4, outfile = "")
    clusterExport(cl, c("selection_bound","X","Y","W","N","sample","j","auglag","grad",
                        "jacobian","results_pop"),
                  envir = environment())
  } else {
    cl <- nnodes
  }

  out <- pbsapply(X=1:5e3, cl = cl, simplify=T, FUN=function(x) {
    # Sample data
    n <- sample[j]
    id <- sample(N, n, replace=F, prob=NULL)
    Ys <- Y[id]; Xs <- X[id,]; Ws <- W[id,]

    # Constraints
    cons <- list(
      list('RESP', 0.07),
      list('COVMEAN', Xs, 0.5)
    )

    # Compute sample confidence interval
    results_sample <- selection_bound(y=Ys, x=Xs, w=Ws, L0l=0.05, L0u=0.15, L1=3)

    # Coverage
    coverage <- (results_sample$ci[1] <= results_pop$interval[1]) &
                (results_sample$ci[2] >= results_pop$interval[2])

    return(coverage)
  })

  if (get_os() == "windows") {
    stopCluster(cl)
  }

  main[j] <- mean(out)
  print(paste("Coverage frequency of", main[j], "for sample size", sample[j]))
}
