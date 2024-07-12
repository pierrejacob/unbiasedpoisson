#'@rdname sample_epave
#'@title Sample EPAVE
#'@description Sample estimator of asymptotic variance in the CLT for MCMC averages.
#'
#'@export
sample_epave <- function(single_kernel, coupled_kernel, rinit, nmcmc, burnin, D, h, state_x_0, max_iterations = Inf){
  state <- rinit()
  # sum of h evaluations
  sumh <- 0
  # sum of h^2 evaluations
  sumh2 <- 0
  nh <- 0 
  # sum of G evaluations
  sumG <- 0
  # sum of h * G evaluations
  sumhG <- 0
  # number of G evaluations
  nG <- 0
  # cost of G evaluations
  costG <- 0
  for (i in 1:nmcmc){
    state <- single_kernel(state)
    if (i > burnin){
      sumh <- sumh + h(state$position)
      sumh2 <- sumh2 + h(state$position)^2
      nh <- nh + 1
      if (i %% D == 0){
        G_ <- sample_unbiasedfishy(coupled_kernel, h, state, state_x_0, max_iterations = max_iterations)
        costG <- costG + G_$meetingtime * 2 # cost = 2 per iteration until meeting time
        sumG <- sumG + G_$estimator
        sumhG <- sumhG + h(state$position) * G_$estimator
        nG <- nG + 1
      }
    }
  }
  ## mean of h under pi
  meanh <- sumh / nh
  ## variance of h under pi:
  vMC <- sumh2 / nh - (meanh)^2
  vEPAVE <- -vMC + 2 * (sumhG / nG - meanh * sumG / nG)
  ## cost of estimator
  costG <- costG
  costEPAVE <- nmcmc + costG
  return(list(vEPAVE = vEPAVE, vMC = vMC, meanh = meanh, costEPAVE = costEPAVE, costG = costG, nG = nG))
}

