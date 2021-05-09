selection_bound <- function(y, x, w, L0, L1, cons=NULL, theta=NULL, alpha=0.05) {

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
  if (is.null(theta)) theta <- rep(0,p+1)

  if (!is.null(cons)) {
    # Select critical values
    k <- length(cons)
    zstat2 <- qnorm(1-alpha/4)
    zstat1 <- qnorm(1-alpha/(4*k))

    # Find an initial feasible theta
    theta <- auglag(x0=theta, fn=qloss, gr=qlossgr, lower=c(log(1/L0), rep(log(1/L1),p)),
                    upper=c(log(L0), rep(log(L1),p)), hin=NULL, hinjac=NULL,
                    localsolver=c("SLSQP"), nl.info=F,
                    control=list("xtol_rel"=1e-8, "ftol_rel"=1e-10, "maxeval"=-1))$par
  } else {
    zstat2 <- qnorm(1-alpha/2)
    hin <- NULL
    hinjac <- NULL
  }

  # Compute bounds
  for (bound in c("max","min")) {
    s <- 2*as.numeric(bound=="min")-1

    results <- auglag(x0=theta, fn=fn, gr=gr, lower=c(log(1/L0), rep(log(1/L1),p)),
                      upper=c(log(L0), rep(log(L1),p)), hin=hin, hinjac=hinjac,
                      localsolver=c("SLSQP"), nl.info=F,
                      control=list("xtol_rel"=1e-8, "ftol_rel"=1e-10, "maxeval"=-1))

    print(results)

    if (any(hin(results$par) < -1e-8)) stop("Solution is not feasible.")
    assign(paste("theta_",bound,sep=""), results$par)
  }

  inv_wgt <- 1+exp(-w%*%theta_min)
  model_min <- summary(lm(y ~ x, weights=inv_wgt))
  ci_lower <- model_min$coefficients[2,1] - zstat2*model_min$coefficients[2,2]

  inv_wgt <- 1+exp(-w%*%theta_max)
  model_max <- summary(lm(y ~ x, weights=inv_wgt))
  ci_upper <- model_max$coefficients[2,1] + zstat2*model_max$coefficients[2,2]

  return(c(ci_lower, ci_upper))
}
