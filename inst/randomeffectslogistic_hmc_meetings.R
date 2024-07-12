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
source("inst/hmc_functions.R")




# sample_meetingtime(HMC, CoupledHMC, rinit, max_iterations = 1000)

# # meeting times
lag <- 500
nrep <- 500
mts <- foreach(irep = 1:nrep, .combine = c) %dorng% {
  res <- sample_meetingtime(HMC, CoupledHMC, rinit, lag = lag)
  res$meetingtime
}

save(mts, nrep, lag, file = "output/randomeffectslogistic.hmc.meetings.RData")
load(file = "output/randomeffectslogistic.hmc.meetings.RData")

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
