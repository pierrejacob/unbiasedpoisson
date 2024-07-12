rm(list = ls())
library(unbiasedpoisson)
setmytheme()
library(doParallel)
library(doRNG)
library(dplyr)
# register parallel cores
registerDoParallel(cores = detectCores()-2)
# set RNG seed
set.seed(1)

## import functions
source("inst/randomeffectslogistic_model.R")
source("inst/mala_functions.R")

### run MALA 
# nmcmc <- 10000
# state <- rinit()
# betahistory <- rep(0, nmcmc)
# betahistory[1] <- state$position[506]
# naccepts <- 0
# for (i in 2:nmcmc){
#   state <- MALA(state)
#   naccepts <- naccepts + state$accept
#   betahistory[i] <- h(state$position)
#   if (i %% 100 == 0){
#     cat("Iteration ", i, "\n")
#     cat("Accept rate:", naccepts/i*100, "%\n")
#   }
# }
# matplot(betahistory, type = "l", lty = 1, col = 1, xlab = "Iteration", ylab = "beta[506]")
# sample_meetingtime(MALA, CoupledMALA, rinit, max_iterations = 10000)

lag <- 10000
nrep <- 500
mts <- foreach(irep = 1:nrep, .combine = c) %dorng% {
  res <- sample_meetingtime(MALA, CoupledMALA, rinit, lag = lag)
  res$meetingtime
}

save(mts, nrep, lag, file = "output/randomeffectslogistic.mala.meetings.RData")
load(file = "output/randomeffectslogistic.mala.meetings.RData")

hist(mts-lag)
mean(mts-lag)


# 
## plot TV upper bounds
niterations <- max(mts)
ubounds <- sapply(1:niterations, function(t) tv_upper_bound(mts, lag, t))
g_tvbounds <- qplot(x = 1:niterations, y = ubounds, geom = "line") +
  ylab("TV distance") + xlab("t") + geom_hline(yintercept = 1)
g_tvbounds <- g_tvbounds + scale_y_log10(breaks = c(0, 0.001, 0.01, 0.1,1)) + geom_rangeframe()
g_tvbounds

