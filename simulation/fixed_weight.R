fixed_weight <- function(y, x, w, theta_min, theta_max, alpha=0.05) {

  # Ensure data are formatted correctly
  y <- as.vector((y - mean(y))/sd(y))
  x <- as.matrix(x); x <- apply(x, 2, function(x) (x - mean(x))/sd(x))
  w <- as.matrix(w); w <- apply(w, 2, function(w) (w - mean(w))/sd(w))

  # Sample size
  n <- length(y)
  if (n != nrow(x) | n != nrow(w)) stop('Rows are unequal.')

  w <- cbind(rep(1,n),w)

  zstat2 <- qnorm(1-alpha/2)

  data <- data.frame(x=x,y=y)

  inv_wgt_min <- 1+exp(-w%*%theta_min)
  model_min <- summary(svyglm(y ~ x, design=svydesign(ids=~0, weights=~inv_wgt_min, data=data)))

  inv_wgt_max <- 1+exp(-w%*%theta_max)
  model_max <- summary(svyglm(y ~ x, design=svydesign(ids=~0, weights=~inv_wgt_max, data=data)))

  out <- list(
    interval = c(model_min$coefficients[2,1], model_max$coefficients[2,1]),
    ci = c(model_min$coefficients[2,1] - zstat2*model_min$coefficients[2,2],
           model_max$coefficients[2,1] + zstat2*model_max$coefficients[2,2])
  )

  return(out)
}
