# Inequality constraints
hin <- function(theta) {
  inv_wgt <- 1+exp(-w%*%theta)

  main <- sapply(X=cons, FUN=function(item) {
    if (item[[1]] == 'RESP') {
      resp <- item[[2]]
      rvar <- var(inv_wgt-1/resp)

      out <- -mean(inv_wgt-1/resp)^2 + zstat^2*rvar/n
    }

    else if (item[[1]] == 'COVMEAN') {
      ccov <- item[[2]]
      cmean <- item[[3]]
      cvar <- var(inv_wgt*(ccov-cmean))

      out <- -mean(inv_wgt*(ccov-cmean))^2 + zstat^2*cvar/n
    }

    else stop(paste(item[[1]], 'is an invalid option.'))

    return(out)
  })

  return(matrix(main, ncol=1))
}

# Jacobian of inequality constraints
hinjac <- function(theta) {
  return(jacobian(hin, theta))
}

# Objective function
fn <- function(theta) {
  inv_wgt <- as.vector(1+exp(-w%*%theta))
  beta <- lm.wfit(x=cbind(rep(1,n),x), y=y, w=inv_wgt)$coefficients[2]
  return(s*beta)
}

# Gradient of objective function
gr <- function(theta) {
  return(grad(fn, theta))
}

# Quadratic loss function for finding feasible region
qloss <- function(theta) {
  return(-sum(hin(theta)-abs(hin(theta))))
}

# Gradient of quadratic loss function for finding feasible region
qlossgr <- function(theta) {
  return(grad(qloss, theta))
}
