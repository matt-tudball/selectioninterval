selection_bound <- function(y, x, w, L0l, L0u, L1, cons=NULL, theta=NULL, alpha=0.05) {

  # Ensure data are formatted correctly
  y <- as.vector(y)
  x <- as.matrix(x)
  w <- as.matrix(w)

  # Sample size
  n <- length(y)
  if (n != nrow(x) | n != nrow(w)) stop('Rows are unequal.')

  # Number of parameters
  p <- ncol(w)
  w <- cbind(rep(1,n),w)

  # Select theta if none provided
  if (is.null(theta)) theta <- c(log(L0l*L0u/((1-L0l)*(1-L0u)))/2, rep(0,p))

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

  if (!is.null(cons)) {
    # Select critical values
    k <- length(cons)
    zstat2 <- qnorm(1-alpha/4)
    zstat1 <- qnorm(1-alpha/(4*k))

    # Inequality constraints
    hin <- function(theta) {
      inv_wgt <- 1+exp(-w%*%theta)

      main <- sapply(X=cons, FUN=function(item) {
        if (item[[1]] == 'RESP') {
          resp <- item[[2]]
          rvar <- var(inv_wgt-1/resp)

          out <- -mean(inv_wgt-1/resp)^2 + zstat1^2*rvar/n
        }

        else if (item[[1]] == 'COVMEAN') {
          ccov <- item[[2]]
          cmean <- item[[3]]
          cvar <- var(inv_wgt*(ccov-cmean))

          out <- -mean(inv_wgt*(ccov-cmean))^2 + zstat1^2*cvar/n
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

    # Quadratic loss function for finding feasible region
    qloss <- function(theta) {
      return(-sum(hin(theta)-abs(hin(theta))))
    }

    # Gradient of quadratic loss function for finding feasible region
    qlossgr <- function(theta) {
      return(grad(qloss, theta))
    }

    # Find an initial feasible theta
    suppressMessages(
    theta <- auglag(x0=theta, fn=qloss, gr=qlossgr, lower=c(log(L0l/(1-L0l)),
                    rep(log(1/L1),p)), upper=c(log(L0u/(1-L0u)), rep(log(L1),p)),
                    hin=NULL, hinjac=NULL, localsolver=c("SLSQP"), nl.info=F,
                    control=list("xtol_rel"=1e-8, "ftol_rel"=1e-10, "maxeval"=-1))$par
    )
  } else {
    zstat2 <- qnorm(1-alpha/2)

    hin <- NULL

    hinjac <- NULL
  }

  # Compute bounds
  for (bound in c("max","min")) {
    s <- 2*as.numeric(bound=="min")-1

    suppressMessages(
      results <- auglag(x0=theta, fn=fn, gr=gr, lower=c(log(L0l/(1-L0l)),
                      rep(log(1/L1),p)), upper=c(log(L0u/(1-L0u)), rep(log(L1),p)),
                      hin=hin, hinjac=hinjac, localsolver=c("MMA"), nl.info=F,
                      control=list("xtol_rel"=1e-8, "ftol_rel"=1e-10, "maxeval"=-1))
    )
    if (!is.null(cons)) {
      if (any(hin(results$par) < -1e-8)) print(hin(results$par)) #stop("Solution is not feasible.")
    }
    assign(paste("theta_",bound,sep=""), results$par)
  }

  inv_wgt <- 1+exp(-w%*%theta_min)
  model_min <- summary(lm(y ~ x, weights=inv_wgt))

  inv_wgt <- 1+exp(-w%*%theta_max)
  model_max <- summary(lm(y ~ x, weights=inv_wgt))

  out <- list(
    interval = c(model_min$coefficients[2,1], model_max$coefficients[2,1]),
    ci = c(model_min$coefficients[2,1] - zstat2*model_min$coefficients[2,2],
           model_max$coefficients[2,1] + zstat2*model_max$coefficients[2,2])
  )

  return(out)
}
