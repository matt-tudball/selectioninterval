iv_selection_bound <- function(y, x, z, w, L0l, L0u, L1, cons=NULL, theta=NULL, alpha=0.05) {

  # Ensure data are formatted correctly
  y <- as.vector(y)
  x <- as.matrix(x)
  z <- as.matrix(z)
  w <- as.matrix(w); w <- apply(w, 2, function(w) (w - mean(w))/sd(w))

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
    beta <- ivreg.fit(x=cbind(1,x), y=y, z=cbind(1,z), weights=inv_wgt)$coefficients[2]
    return(s*beta)
  }

  # Gradient of objective function
  gr <- function(theta) {
    return(grad(fn, theta))
  }

  # Set lower and upper bounds
  lower <- c(log(L0l/(1-L0l)), rep(log(1/L1),p))
  upper <- c(log(L0u/(1-L0u)), rep(log(L1),p))

  for (item in cons) {
    if (item[[1]] == "DIREC") {
      if (item[[3]] == "+") lower[(item[[2]]+1)] <- 0
      else if (item[[3]] == "-") upper[(item[[2]]+1)] <- 0
      else stop(paste(item[[3]], "is an invalid direction."))
    }
  }

  cons <- cons[sapply(X=cons, FUN=function(item) {
    return(item[[1]] %in% c("RESP","COVMEAN"))
  })]

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
          ccov <- w[,item[[2]]+1]
          cmean <- item[[3]]
          cvar <- var(inv_wgt*(ccov-cmean))
          out <- -mean(inv_wgt*(ccov-cmean))^2 + zstat1^2*cvar/n
        }

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
      theta <- auglag(x0=theta, fn=qloss, gr=qlossgr, lower=lower, upper=upper,
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
      global_solve <- nloptr(x0=theta, eval_f=fn, lb=lower, ub=upper,
                             eval_g_ineq=hin, opts=list("xtol_rel"=1e-2, "maxeval"=1e5,
                             "algorithm"="NLOPT_GN_ISRES","print_level"=1))
    )

    theta <- global_solve$solution

    suppressMessages(
      results <- auglag(x0=theta, fn=fn, gr=gr, lower=lower, upper=upper,
                        hin=hin, hinjac=hinjac, localsolver=c("SLSQP"), nl.info=F,
                        control=list("xtol_rel"=1e-8, "ftol_rel"=1e-10, "maxeval"=-1))
    )
    if (!is.null(cons)) {
      if (any(hin(results$par) < -1e-8)) stop("Solution is not feasible.")
    }
    assign(paste("theta_",bound,sep=""), results$par)
  }

  data <- data.frame(z=z,x=x,y=y)
  xnames <- colnames(data)[startsWith(colnames(data),'x')]
  znames <- colnames(data)[startsWith(colnames(data),'z')]

  formula <- as.formula(paste('y ~', paste(xnames, collapse='+'), '|', paste(znames, collapse='+'),sep=''))

  inv_wgt_min <- as.vector(1+exp(-w%*%theta_min))
  model_min <- summary(ivreg(formula=formula, data=data, weights=inv_wgt_min))

  inv_wgt_max <- as.vector(1+exp(-w%*%theta_max))
  model_max <-  summary(ivreg(formula=formula, data=data, weights=inv_wgt_max))

  out <- list(
    theta_min = theta_min,
    theta_max = theta_max,
    interval = c(model_min$coefficients[2,1], model_max$coefficients[2,1]),
    ci = c(model_min$coefficients[2,1] - zstat2*model_min$coefficients[2,2],
           model_max$coefficients[2,1] + zstat2*model_max$coefficients[2,2])
  )

  return(out)
}
