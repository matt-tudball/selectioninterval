#' @name selection_bound
#' @title Sensitivity analysis for selection bias
#' @description Given a set of sensitivity parameters and constraints, computes an
#' upper and lower bound for an inverse probability weighted regression estimate
#'
#' @param y Outcome (vector)
#' @param x Explanatory variables (matrix)
#' @param w Selection variables (matrix)
#' @param z Optional instrumental variables (matrix)
#' @param L0l Lower bound for the probability of sample selection for the average observation
#' @param L0u Upper bound for the probability of sample selection for the average observation
#' @param L1 Odds ratio of sample selection for a one unit (binary) or one standard
#' deviation (continuous) increase in a variable.
#' @param cons List of constraints to be applied: \code{RESP} is a response rate
#' constraint; \code{COVMEAN} is a covariate mean constraint; \code{DIREC} is a
#' directionality constraint.
#' @param theta Optional starting parameter for the global optimiser.
#' @param alpha Significance level for confidence interval.
#' @param opts Optional list of options for the global optimiser.
#'
#' @return Named list of objects: \code{theta_min} (\code{theta_max}) is the parameters
#' for the lower (upper) bound; \code{interval} is a vector containing the lower and upper
#' bounds; \code{ci} is a vector containing the \code{alpha}-level confidence interval.
#'
#' @importFrom stats lm.wfit
#' @importFrom survey svydesign svyglm
#' @importFrom nloptr nloptr auglag
#' @importFrom numDeriv grad jacobian
#' @importFrom AER ivreg ivreg.fit
#'
#' @examples
#'

selection_bound <- function(y, x, w, z=NULL, L0l, L0u, L1, cons=NULL, theta=NULL, alpha=0.05, opts=NULL) {

  # Ensure data are formatted correctly
  y <- as.vector(y)
  x <- as.matrix(x)
  if(!is.null(z)) z <- as.matrix(z)
  w <- as.matrix(w); w <- apply(w, 2, function(v) (v - mean(v))/sd(v))

  # Sample size
  n <- length(y)
  if (n != nrow(x) | n != nrow(w)) stop('Rows are unequal.')

  # Number of parameters
  p <- ncol(w)
  w <- cbind(rep(1,n),w)

  # Replace defaults with user-supplied options
  global_opts <- list(algorithm="NLOPT_GN_ISRES", maxeval=5e4, xtol_rel=1e-2, print_level=0)
  if (!is.null(opts)) {
    for (j in 1:length(opts)) {
      if(names(opts[j]) == "algorithm") global_opts$algorithm <- opts[[j]]
      else if(names(opts[j]) == "maxeval") global_opts$maxeval <- opts[[j]]
      else if(names(opts[j]) == "xtol_rel") global_opts$xtol_rel <- opts[[j]]
      else if(names(opts[j]) == "print_level") global_opts$print_level <- opts[[j]]
      else stop(paste(names(opts[j]), "is not a valid option."))
    }
  }

  # Select theta if none provided
  if (is.null(theta)) theta <- c(log(L0l*L0u/((1-L0l)*(1-L0u)))/2, rep(0,p))

  # Objective function (regression or instrumental variable)
  if(!is.null(z)) {
    fn <- function(theta) {
      inv_wgt <- as.vector(1+exp(-w%*%theta))
      beta <- ivreg.fit(x=cbind(1,x), y=y, z=cbind(1,z), weights=inv_wgt)$coefficients[2]
      return(s*beta)
    }
  } else {
    fn <- function(theta) {
      inv_wgt <- as.vector(1+exp(-w%*%theta))
      beta <- lm.wfit(x=cbind(rep(1,n),x), y=y, w=inv_wgt)$coefficients[2]
      return(s*beta)
    }
  }

  # Numerical radient of objective function
  gr <- function(theta) {
    return(grad(fn, theta))
  }

  # Set lower and upper bounds
  lower <- c(log(L0l/(1-L0l)), rep(log(1/L1),p))
  upper <- c(log(L0u/(1-L0u)), rep(log(L1),p))

  # Modify the bounds based on DIREC constraints
  for (item in cons) {
    if (item[[1]] == "DIREC") {
      if (item[[3]] == "+") lower[(item[[2]]+1)] <- 0
      else if (item[[3]] == "-") upper[(item[[2]]+1)] <- 0
      else stop(paste(item[[3]], "is an invalid direction."))
    }
  }

  # Keep RESP and COVMEAN constraints
  ind <- sapply(X=cons, FUN=function(item) { return(item[[1]] %in% c("RESP","COVMEAN")) })
  if (any(ind)) { cons <- cons <- cons[ind] } else { cons <- NULL }

  # Select critical values and constraint functions (if needed)
  if (!is.null(cons)) {
    # Select critical values
    k <- length(cons)
    zstat2 <- qnorm(1-alpha/4)
    zstat1 <- qnorm(1-alpha/(4*k))

    # Inequality constraints
    hin <- function(theta) {
      inv_wgt <- as.vector(1+exp(-w%*%theta))

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
        
        else if (item[[1]] == 'NEGCON') {
          # Variable one
          w1 <- item[[2]]
          # Variable two
          w2 <- item[[3]]
          # Checking correlation (want it to be zero)
          ncorr <- lm.wfit(x=w1,y=w2,w=inv_wgt)$coefficients[2]
          # Variance of this correlation
          nvar <- mean(lm.wfit(x=w1,y=w2,w=inv_wgt)$residuals^2)/sum(w2^2)
          # Checking whether constraints hold
          out <- -ncorr^2 + zstat1^2*nvar
        }

        else if (item[[1]] == 'NEGCON') {
          w1 <- item[[2]]
          w2 <- item[[3]]
          nx <- as.matrix(cbind(1,w1))
          nwgt <- diag(inv_wgt)
          ninv <- solve(t(nx)%*%nwgt%*%nx)

          ncoef <- ninv%*%t(nx)%*%nwgt%*%w2

          nres <- as.vector((w2 - nx%*%ncoef)^2)
          nvar <- ninv%*%t(nx)%*%nwgt%*%diag(nres)%*%t(nwgt)%*%nx%*%ninv

          # Checking whether constraints hold
          out <- ncoef[2]^2 + zstat1^2*nvar[2,2]
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
      theta0 <- auglag(x0=theta, fn=qloss, gr=qlossgr, lower=lower, upper=upper,
                      hin=NULL, hinjac=NULL, localsolver=c("SLSQP"), nl.info=F,
                      control=list("xtol_rel"=1e-8, "ftol_rel"=1e-10, "maxeval"=-1))$par
    )

    if (qloss(theta0) > 0) stop(
      "Could not find a feasible starting value. Try increasing the sensitivity parameters
      or including more variables in the weight model."
    )
  } else {
    theta0 <- theta

    zstat2 <- qnorm(1-alpha/2)

    hin <- NULL

    hinjac <- NULL
  }

  # Compute bounds
  for (bound in c("max","min")) {
    s <- 2*as.numeric(bound=="min")-1

    # Global optimiser to find approximate solution
    suppressMessages(
      global_solve <- nloptr(x0=theta0, eval_f=fn, lb=lower, ub=upper, eval_g_ineq=hin,
                             opts=global_opts)
    )

    theta1 <- global_solve$solution

    # Local optimiser to refine approximate solution
    suppressMessages(
      results <- auglag(x0=theta1, fn=fn, gr=gr, lower=lower, upper=upper,
                        hin=hin, hinjac=hinjac, localsolver=c("SLSQP"), nl.info=F,
                        control=list("xtol_rel"=1e-8, "ftol_rel"=1e-10, "maxeval"=-1))
    )

    # Check if final solution is feasible
    if (!is.null(cons)) {
      if (any(hin(results$par) < -1e-8)) stop(
        "Solution is not feasible. This error originates with the nloptr package.
        Running the function again sometimes resolves this error."
      )
    }
    assign(paste("theta_",bound,sep=""), results$par)
  }

  # Compute main objects
  if (!is.null(z)) {
    data <- data.frame(z=z,x=x,y=y)
    xnames <- colnames(data)[startsWith(colnames(data),'x')]
    znames <- colnames(data)[startsWith(colnames(data),'z')]
    formula <- as.formula(paste('y ~', paste(xnames, collapse='+'), '|', paste(znames, collapse='+'),sep=''))

    inv_wgt_min <- as.vector(1+exp(-w%*%theta_min))
    model_min <- summary(ivreg(formula=formula, data=data, weights=inv_wgt_min))

    inv_wgt_max <- as.vector(1+exp(-w%*%theta_max))
    model_max <-  summary(ivreg(formula=formula, data=data, weights=inv_wgt_max))
  } else {
    data <- data.frame(x=x,y=y)

    inv_wgt_min <- 1+exp(-w%*%theta_min)
    model_min <- summary(svyglm(y ~ x, design=svydesign(ids=~0, weights=inv_wgt_min, data=data)))

    inv_wgt_max <- 1+exp(-w%*%theta_max)
    model_max <- summary(svyglm(y ~ x, design=svydesign(ids=~0, weights=inv_wgt_max, data=data)))
  }

  # Return optimal solutions, interval and confidence interval
  out <- list(
    theta_min = theta_min,
    theta_max = theta_max,
    interval = c(model_min$coefficients[2,1], model_max$coefficients[2,1]),
    ci = c(model_min$coefficients[2,1] - zstat2*model_min$coefficients[2,2],
           model_max$coefficients[2,1] + zstat2*model_max$coefficients[2,2])
  )

  return(out)
}
