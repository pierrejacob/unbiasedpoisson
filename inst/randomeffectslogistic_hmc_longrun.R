rm(list = ls())
library(unbiasedpoisson)
library(doParallel)
library(doRNG)
set.seed(1)
setmytheme()
registerDoParallel(cores = detectCores() - 2)
source("inst/randomeffectslogistic_model.R")
source("inst/hmc_functions.R")


nmcmc <- 100000
nchains <- 100
# filepath <- "~/Dropbox/UnbiasedPoissonNumerics"
filepath <- ""


history <- foreach(ichain = 1:nchains, .combine = cbind) %dorng% {
  history_onechain <- rep(0, nmcmc)
  state <- rinit()
  for (imcmc in 1:nmcmc){
    state <- HMC(state)
    history_onechain[imcmc] <- state$position[506]
  }
  history_onechain
}

history.df <- reshape2::melt(history)
names(history.df) <- c("iteration", "chain", "value")
history.df$chain <- rep(1:nchains, each = nmcmc)

save(nmcmc, nchains, stepsize, Lmax, history.df, file = file.path(filepath, "randomeffects.hmc.longrun.RData"))
load(file = file.path(filepath, "randomeffects.hmc.longrun.RData"))

tail(history.df)

library(tidyverse)
ggplot(history.df %>% filter(chain <= 10), aes(x = iteration, y = value, color = chain)) + geom_line() + theme_minimal() +
  theme(legend.position = "none")
# 
# ## run HMC
# state <- rinit()
# 
# betahistory <- rep(0, nmcmc)
# betahistory[1] <- state$position[506]
# naccepts <- 0
# for (i in 2:nmcmc){
#   state <- HMC(state)
#   naccepts <- naccepts + state$accept
#   betahistory[i] <- state$position[506]
#   if (i %% 100 == 0){
#     cat("Iteration ", i, "\n")
#     cat("Accept rate:", naccepts/i*100, "%\n")
#   }
# }
# matplot(betahistory, type = "l", lty = 1, col = 1, xlab = "Iteration", ylab = "beta[506]")
