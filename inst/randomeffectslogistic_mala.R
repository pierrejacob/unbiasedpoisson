rm(list = ls())
library(unbiasedpoisson)
library(doParallel)
library(doRNG)
set.seed(1)
setmytheme()
registerDoParallel(cores = detectCores() - 2)
source("inst/randomeffectslogistic_model.R")
# plot(invmassdiag, diag(postcov))

source("inst/mala_functions.R")


## define diagonal proposal covariance
Sigma_prop <- diag(postcov)
Sigma_prop_chol <- sqrt(Sigma_prop)
# inv_Sigma_prop_chol <- 1/Sigma_prop_chol
## define scalar stepsize
stepsize <- target_dim^{-1/6}

stepsize <- 0.2

## run MALA
state <- rinit()
nmcmc <- 10000
betahistory <- rep(0, nmcmc)
betahistory[1] <- state$position[506]
naccepts <- 0
for (i in 2:nmcmc){
  state <- MALA(state)
  naccepts <- naccepts + state$accept
  betahistory[i] <- state$position[506]
  if (i %% 100 == 0){
    cat("Iteration ", i, "\n")
    cat("Accept rate:", naccepts/i*100, "%\n")
  }
}

matplot(betahistory, type = "l", lty = 1, col = 1, xlab = "Iteration", ylab = "beta[506]")

# ## run coupled MALA
# lag <- 5e3
# nrep <- 100
# ccs <- foreach(irep = 1:nrep) %dorng% {
#   sample_coupled_chains(MALA, CoupledMALA, rinit, lag = lag)
# }
# 
# ##
# mts <- sapply(X = ccs, function(x) x$meetingtime, simplify = T)
# summary(mts)
# plot(ccs[[which.max(mts)]]$samples1[,506], type = "l", lty = 1, col = 1, xlab = "Iteration", ylab = "beta[506]")
# lines(x = (lag+1):(dim(ccs[[which.max(mts)]]$samples2)[1]+lag), y = ccs[[which.max(mts)]]$samples2[1:dim(ccs[[which.max(mts)]]$samples2)[1],506], col = 'red')


# 
# niterations <- max(mts)
# ubounds <- sapply(1:niterations, function(t) tv_upper_bound(mts, lag, t))
# g_tvbounds <- qplot(x = 1:niterations, y = ubounds, geom = "line") +
#   ylab("TV distance") + xlab("t") + geom_hline(yintercept = 1)
# g_tvbounds <- g_tvbounds + scale_y_log10(breaks = c(0, 0.001, 0.01, 0.1,1)) + geom_rangeframe() 
# g_tvbounds
# 
# 
