rm(list=ls())
require(ggplot2); require(parallel); require(devtools); require(pbapply)
setwd("C:/Users/ow18301/OneDrive - University of Bristol/MyFiles-Migrated/Documents")

# Load functions
load_all(path='IEASBOS_FILES_/CODE/selectioninterval')

# Seed
set.seed(2021)

# Population size
N <- 1e6

# Population data
#X <- as.matrix(rbeta(N,2,5))
X <- as.matrix(rnorm(N,0,1))
#X <- as.matrix(rexp(N,1))
Y <- as.vector(rnorm(N,0,1))
W <- as.matrix(cbind(X,Y))

mycons <- list(list('DIREC',1,'+'), list('RESP',0.15))

# Population confidence interval
results_pop <- selection_bound(y=Y, x=X, w=W, L0l=0.1, L0u=0.2, L1=3, cons=mycons)

# Simulation
sample <- c(seq(10,50,5), seq(100,1000,100))#, seq(2000,10000,1000), seq(12000,20000,2000))
main <- rep(NA, length(sample))

for (j in 1:length(sample)) {
  if (get_os() == "windows") {
    cl <- makeCluster(4, outfile = "")
    clusterExport(cl, c("selection_bound","X","Y","W","N","sample","j","auglag","grad",
                        "jacobian","results_pop","nloptr","svydesign","svyglm","mycons"),
                  envir = environment())
  } else {
    cl <- nnodes
  }

  out <- t(pbsapply(X=1:5e3, cl = cl, simplify=T, FUN=function(x) {
    # Sample data
    n <- sample[j]
    id <- sample(N, n, replace=T, prob=NULL)
    Ys <- Y[id]; Xs <- X[id,]; Ws <- W[id,]

    # Compute sample confidence interval
    # In rare instances the underlying optimiser (nlopt) will fail to find a solution.
    # As far as I can tell, this is a computational not a statistical problem.
    results_sample <- try(selection_bound(y=Ys, x=Xs, w=Ws, L0l=0.1, L0u=0.2, L1=3, cons=mycons))
    if (inherits(results_sample, "try-error")) {
      coverage <- NA
    } else { coverage <- (results_sample$ci[1] <= results_pop$interval[1]) &
      (results_sample$ci[2] >= results_pop$interval[2])
    }

    return(coverage)
  }))

  if (get_os() == "windows") {
    stopCluster(cl)
  }

  main[j] <- mean(out, na.rm=T)
  print(paste("Coverage is ", main[j]," for sample size ", sample[j], sep=""))
}

saveRDS(main, file="IEASBOS_FILES_/SIMULATIONS/DATA/coverage_scenario3.rds")

# Make the plot
data <- data.frame(x=sample, y=main)
plot <- ggplot(data=data, aes(x=x, y=y)) +
  geom_line() +
  geom_point() +
  xlab("Sample size") +
  ylab("Coverage frequency") +
  scale_y_continuous(breaks=seq(0,1,0.2)) +
  geom_hline(yintercept=0.95, color="grey42")
#ggsave(filename="IEASBOS_FILES_/SIMULATIONS/FIGURES/coverage.pdf",plot=plot)

