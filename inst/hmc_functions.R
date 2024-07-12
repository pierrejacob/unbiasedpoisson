source("inst/mala_functions.R")

## HMC with
## - random stepsize, drawn from Exponential(epsilon)
## - random number of leapfrog steps L ~ Uniform(1, Lmax)
## - diagonal mass matrix M = Sigma_prop, and Sigma_prop_chol = sqrt(Sigma_prop)

## The function requires
## - target: a function that takes a vector as input and returns the log-density of the target
## - gradtarget: a function that takes a vector as input and returns the gradient of the log-density of the target
## - target_dim: an integer specifying the dimension of the target
## - Sigma_prop: a vector of length target_dim, specifying the diagonal of the mass matrix
## - Sigma_prop_chol: a vector of length target_dim, specifying the square root of the diagonal of the mass matrix
## - Lmax: an integer specifying the maximum number of leapfrog steps
## - stepsize: a positive scalar specifying the stepsize

## The argument 'state' must be a list with
## - position: a vector of length target_dim, specifying the current position of the chain
## - current_pdf: a scalar, specifying the current log-density of the target at the current position
## - grad_pdf: a vector of length target_dim, specifying the gradient of the log-density of the target at the current position

HMC <- function(state) {
  nleapfrogsteps <- sample.int(Lmax, 1)
  if (nleapfrogsteps == 1){
    return(MALA(state))
  } else {
    eps <- rexp(rate = 1/stepsize, n = 1)
    # draw momentum
    ## mass matrix M is inverse of Sigma_prop
    initial_momentum <- rnorm(target_dim, 0, sd=1/Sigma_prop_chol)
    position <- state$position
    g <- state$grad_pdf
    # numerical error
    num_error <- FALSE
    # leap frog integrator
    momentum <- initial_momentum + eps * g / 2
    for (step in 1:nleapfrogsteps){
      position <- position + eps * momentum * Sigma_prop
      g <- gradtarget(position)
      if (any(is.na(g)) || any(is.infinite(g))){
        num_error <- TRUE
        break
      }
      if (step != nleapfrogsteps){
        momentum <- momentum + eps * g
      }
    }
    if (num_error){
      return(list(position = state$position, current_pdf = state$current_pdf, grad_pdf = state$grad_pdf, accept = FALSE))
    }
    momentum <- momentum + eps * g / 2
    proposed_pdf <- target(position)
    accept_ratio <- proposed_pdf - state$current_pdf
    # the acceptance ratio also features the "kinetic energy" term of the extended target
    accept_ratio <- accept_ratio + (-0.5 * sum(momentum^2 * Sigma_prop)) - 
      (-0.5 * sum(initial_momentum^2 * Sigma_prop))
    accept <- FALSE
    logu <- log(runif(1))
    if (is.finite(accept_ratio)){
      accept <- (logu < accept_ratio)
    }
    if (accept){
      return(list(position = position, current_pdf = proposed_pdf, grad_pdf = g, accept = accept))
    } else {
      return(list(position = state$position, current_pdf = state$current_pdf, grad_pdf = state$grad_pdf, accept = accept))
    }
  }
}

## Coupled HMC where, 
## - if L = 1 (MALA), we use a coupled MALA kernel with reflection-maximal coupling of proposal
## - if L > 1, we use common random numbers

## Requirements and arguments are as for the function "HMC"

CoupledHMC <- function(state1, state2) {
  nleapfrogsteps <- sample.int(Lmax, 1)
  if (nleapfrogsteps == 1){
    return(CoupledMALA(state1, state2))
  } else {
    eps <- rexp(rate = 1/stepsize, n = 1)
    # draw momentum
    ## mass matrix M is inverse of Sigma_prop
    initial_momentum <- rnorm(target_dim, 0, sd=1/Sigma_prop_chol)
    position1 <- state1$position
    g1 <- state1$grad_pdf
    position2 <- state2$position
    g2 <- state2$grad_pdf
    # numerical error
    num_error1 <- FALSE
    num_error2 <- FALSE
    # leap frog integrator
    momentum1 <- initial_momentum + eps * g1 / 2
    momentum2 <- initial_momentum + eps * g2 / 2
    for (step in 1:nleapfrogsteps){
      if (!num_error1){
        position1 <- position1 + eps * momentum1 * Sigma_prop
        g1 <- gradtarget(position1)
        num_error1 <- any(is.na(g1)) || any(is.infinite(g1))
      }
      if (!num_error2){
        position2 <- position2 + eps * momentum2 * Sigma_prop
        g2 <- gradtarget(position2)
        num_error2 <- any(is.na(g2)) || any(is.infinite(g2))
      }
      if (step != nleapfrogsteps){
        momentum1 <- momentum1 + eps * g1
        momentum2 <- momentum2 + eps * g2
      }
    }
    momentum1 <- momentum1 + eps * g1 / 2
    momentum2 <- momentum2 + eps * g2 / 2
    logu <- log(runif(1))
    ## accept/reject mechanism for state 1
    accept1 <- FALSE
    if (num_error1){
      newstate1 <- state1
      newstate1$accept <- accept1
    } else {
      proposed_pdf1 <- target(position1)
      accept_ratio1 <- proposed_pdf1 - state1$current_pdf
      accept_ratio1 <- accept_ratio1 + (-0.5 * sum(momentum1^2 * Sigma_prop)) - 
        (-0.5 * sum(initial_momentum^2 * Sigma_prop))
      if (is.finite(accept_ratio1)){
        accept1 <- (logu < accept_ratio1)
      }
      if (accept1){
        newstate1 <- list(position = position1, current_pdf = proposed_pdf1, grad_pdf = g1, accept = accept1)
      } else {
        newstate1 <- state1
        newstate1$accept <- accept1
      }
    }
    ## accept/reject mechanism for state 2
    accept2 <- FALSE
    if (num_error2){
      newstate2 <- state2
      newstate2$accept <- accept2
    } else {
      proposed_pdf2 <- target(position2)
      accept_ratio2 <- proposed_pdf2 - state2$current_pdf
      accept_ratio2 <- accept_ratio2 + (-0.5 * sum(momentum2^2 * Sigma_prop)) - 
        (-0.5 * sum(initial_momentum^2 * Sigma_prop))
      if (is.finite(accept_ratio2)){
        accept2 <- (logu < accept_ratio2)
      }
      if (accept2){
        newstate2 <- list(position = position2, current_pdf = proposed_pdf2, grad_pdf = g2, accept = accept2)
      } else {
        newstate2 <- state2
        newstate2$accept <- accept2
      }
    }
    return(list(state1 = newstate1, state2 = newstate2, 
                identical = FALSE))
  }
}


