# MALA
# X* ~ Normal(mean=X + (stepsize)^2/2 * Sigma_prop * grad_pdf(X), var=stepsize^2 * Sigma_prop)
# in this parametrization, stepsize should scale as dimension^{-1/6}
# requires target, gradtarget, stepsize, Sigma_prop, Sigma_prop_chol, inv_Sigma_prop_chol
# e.g.
# stepsize <- dimension^{-1/6}
# Sigma_prop <- diag(1, dimension, dimension)
# Sigma_prop_chol <- chol(Sigma_prop)
# inv_Sigma_prop_chol <- solve(Sigma_prop_chol)


##  
## The MALA function requires
## - target: a function that takes a vector of length target_dim and returns the log-density of the target at that point
## - gradtarget: a function that takes a vector of length target_dim and returns the gradient of the log-density of the target at that point
## - stepsize: a scalar, specifying the stepsize of the MALA proposal
## - target_dim: an integer specifying the dimension of the target (actually not used in the current function)
## - Sigma_prop: a vector of length target_dim, specifying the diagonal mass matrix
## - Sigma_prop_chol: a vector of length target_dim, specifying the square root of the diagonal mass matrix
## - inv_Sigma_prop_chol: a vector of length target_dim, specifying the inverse of the square root of the diagonal mass matrix

## The argument 'state' must be a list with
## - position: a vector of length target_dim, specifying the current position of the chain
## - current_pdf: a scalar, specifying the current log-density of the target at the current position
## - grad_pdf: a vector of length target_dim, specifying the gradient of the log-density of the target at the current position

## Below you can find a function MALA_fullcov that works with a full covariance for the proposal

MALA <- function(state) {
  eps <- rexp(rate = 1/stepsize, n = 1)
  x <- state$position
  g <- state$grad_pdf
  meanprop <- x + (eps)^2/2 * g * Sigma_prop
  xprop <- rnorm(length(meanprop), meanprop, eps * Sigma_prop_chol)
  prop_pdf <- target(xprop)
  gprop <- gradtarget(xprop)
  if (is.finite(prop_pdf) && all(is.finite(gprop))){
    acceptratio <- prop_pdf - state$current_pdf
    ## see Proposition 1 in Titsias, Optimal Preconditioning and Fisher Adaptive Langevin Sampling
    acceptratio <- acceptratio + 0.5 * sum((x - xprop - (eps^2/4) * gprop * Sigma_prop) * gprop) - 
      0.5 * sum((xprop - x - (eps^2/4) * g * Sigma_prop) * g) 
    if (is.nan(acceptratio)){
      acceptratio <- -Inf
    }
  } else {
    acceptratio <- -Inf
  }
  logu <- log(runif(1))
  accept <- (logu < acceptratio)
  if (accept){
    return(list(position = xprop, current_pdf = prop_pdf, grad_pdf = gprop, accept = accept))
  } else {
    return(list(position = x, current_pdf = state$current_pdf, grad_pdf = state$grad_pdf, accept = accept))
  }
}


# Coupled MALA 
# requires target, gradtarget, stepsize, Sigma_prop, Sigma_prop_chol, inv_Sigma_prop_chol
# (same as MALA function above)
# coupling strategy: reflection-maximal coupling at every step, and common uniform for acceptance
CoupledMALA <- function(state1, state2){
  eps <- rexp(rate = 1/stepsize, n = 1)
  x1 <- state1$position
  g1 <- state1$grad_pdf
  meanprop1 <- x1 + (eps)^2/2 * g1 * Sigma_prop
  x2 <- state2$position
  g2 <- state2$grad_pdf
  meanprop2 <- x2 + (eps)^2/2 * g2 * Sigma_prop
  ## coupled proposals from reflection-maximal coupling
  coupledprops <- rmvnorm_reflectionmax_diag(meanprop1, meanprop2, eps * Sigma_prop_chol)
  xprop1 <- coupledprops$xy[,1]
  prop_pdf1 <- target(xprop1)
  gprop1 <- gradtarget(xprop1)
  if (is.finite(prop_pdf1) && all(is.finite(gprop1))){
    acceptratio1 <- prop_pdf1 - state1$current_pdf
    ## see Proposition 1 in Titsias, Optimal Preconditioning and Fisher Adaptive Langevin Sampling
    acceptratio1 <- acceptratio1 + 0.5 * sum((x1 - xprop1 - (eps^2/4) * gprop1 * Sigma_prop) * gprop1) -
      0.5 * sum((xprop1 - x1 - (eps^2/4) * g1 * Sigma_prop) * g1)
    if (is.nan(acceptratio1)){
      acceptratio1 <- -Inf
    }
  } else {
    acceptratio1 <- -Inf
  }
  xprop2 <- coupledprops$xy[,2]
  prop_pdf2 <- target(xprop2)
  gprop2 <- gradtarget(xprop2)
  if (is.finite(prop_pdf2) && all(is.finite(gprop2))){
    acceptratio2 <- prop_pdf2 - state2$current_pdf
    ## see Proposition 1 in Titsias, Optimal Preconditioning and Fisher Adaptive Langevin Sampling
    acceptratio2 <- acceptratio2 + 0.5 * sum((x2 - xprop2 - (eps^2/4) * gprop2 * Sigma_prop) * gprop2) -
      0.5 * sum((xprop2 - x2 - (eps^2/4) * g2 * Sigma_prop) * g2)
    if (is.nan(acceptratio2)){
      acceptratio2 <- -Inf
    }
  } else {
    acceptratio2 <- -Inf
  }
  logu <- log(runif(1))
  accept1 <- (logu < acceptratio1)
  accept2 <- (logu < acceptratio2)
  if (accept1){
    state1 <- list(position = xprop1, current_pdf = prop_pdf1, grad_pdf = gprop1)
  }
  if (accept2){
    state2 <- list(position = xprop2, current_pdf = prop_pdf2, grad_pdf = gprop2)
  }
  return(list(
    state1 = state1,
    state2 = state2,
    identical = (coupledprops$identical && accept1 && accept2)
  ))
}


## here's the code for MALA with full covariance, 
## assuming that 'x' and 'g' are *row* vectors (1xd matrices), not arrays or column vectors
MALA_fullcov <- function(state) {
  eps <- rexp(rate = 1/stepsize, n = 1)
  x <- state$position
  g <- state$grad_pdf
  meanprop <- x + (eps)^2/2 * g %*% Sigma_prop
  xprop <- fast_rmvnorm_chol(1, meanprop, eps * Sigma_prop_chol)
  prop_pdf <- target(xprop)
  gprop <- gradtarget(xprop)
  if (is.finite(prop_pdf) && all(is.finite(gprop))){
    acceptratio <- prop_pdf - state$current_pdf
    ## see Proposition 1 in Titsias, Optimal Preconditioning and Fisher Adaptive Langevin Sampling
    acceptratio <- acceptratio + 0.5 * (x - xprop - (eps^2/4) * gprop %*% Sigma_prop) %*% t(gprop) - 
      0.5 * (xprop - x - (eps^2/4) * g %*% Sigma_prop) %*% t(g) 
    if (is.nan(acceptratio)){
      acceptratio <- -Inf
    }
  } else {
    acceptratio <- -Inf
  }
  logu <- log(runif(1))
  accept <- (logu < acceptratio)
  if (accept){
    return(list(position = xprop, current_pdf = prop_pdf, grad_pdf = gprop, accept = accept))
  } else {
    return(list(position = x, current_pdf = state$current_pdf, grad_pdf = state$grad_pdf, accept = accept))
  }
}

# Coupled MALA 
# still assuming that 'x' and 'g' are *row* vectors (1xd matrices), not arrays or column vectors
# requires target, gradtarget, stepsize, Sigma_prop, Sigma_prop_chol, inv_Sigma_prop_chol
# (same as MALA function above)
# coupling strategy: reflection-maximal coupling at every step, and common uniform for acceptance
CoupledMALA_fullcov <- function(state1, state2){
  eps <- rexp(rate = 1/stepsize, n = 1)
  x1 <- state1$position
  g1 <- state1$grad_pdf
  meanprop1 <- x1 + (eps)^2/2 * g1 %*% Sigma_prop
  x2 <- state2$position
  g2 <- state2$grad_pdf
  meanprop2 <- x2 + (eps)^2/2 * g2 %*% Sigma_prop
  ## coupled proposals from reflection-maximal coupling
  coupledprops <- rmvnorm_reflectionmax(meanprop1, meanprop2, eps * Sigma_prop_chol, inv_Sigma_prop_chol / eps)
  xprop1 <- matrix(coupledprops$xy[,1], nrow = 1)
  prop_pdf1 <- target(xprop1)
  gprop1 <- gradtarget(xprop1)
  if (is.finite(prop_pdf1) && all(is.finite(gprop1))){
    acceptratio1 <- prop_pdf1 - state1$current_pdf
    ## see Proposition 1 in Titsias, Optimal Preconditioning and Fisher Adaptive Langevin Sampling
    acceptratio1 <- acceptratio1 + 0.5 * (x1 - xprop1 - (eps^2/4) * gprop1 %*% Sigma_prop) %*% t(gprop1) - 
      0.5 * (xprop1 - x1 - (eps^2/4) * g1 %*% Sigma_prop) %*% t(g1) 
    if (is.nan(acceptratio1)){
      acceptratio1 <- -Inf
    }
  } else {
    acceptratio1 <- -Inf
  }
  xprop2 <- matrix(coupledprops$xy[,2], nrow = 1)
  prop_pdf2 <- target(xprop2)
  gprop2 <- gradtarget(xprop2)
  if (is.finite(prop_pdf2) && all(is.finite(gprop2))){
    acceptratio2 <- prop_pdf2 - state2$current_pdf
    ## see Proposition 1 in Titsias, Optimal Preconditioning and Fisher Adaptive Langevin Sampling
    acceptratio2 <- acceptratio2 + 0.5 * (x2 - xprop2 - (eps^2/4) * gprop2 %*% Sigma_prop) %*% t(gprop2) - 
      0.5 * (xprop2 - x2 - (eps^2/4) * g2 %*% Sigma_prop) %*% t(g2) 
    if (is.nan(acceptratio2)){
      acceptratio2 <- -Inf
    }
  } else {
    acceptratio2 <- -Inf
  }
  logu <- log(runif(1))
  accept1 <- (logu < acceptratio1)
  accept2 <- (logu < acceptratio2)
  if (accept1){
    state1 <- list(position = xprop1, current_pdf = prop_pdf1, grad_pdf = gprop1)
  }
  if (accept2){
    state2 <- list(position = xprop2, current_pdf = prop_pdf2, grad_pdf = gprop2)
  }
  return(list(
    state1 = state1,
    state2 = state2,
    identical = (coupledprops$identical && accept1 && accept2)
  ))
}

