rm(list = ls())
library(unbiasedpoisson)
library(doParallel)
library(doRNG)
set.seed(1)
setmytheme()
registerDoParallel(cores = detectCores() - 2)
source("inst/randomeffectslogistic_model.R")
source("inst/mala_functions.R")

# filepath <- "~/Dropbox/UnbiasedPoissonNumerics"
filepath <- ""

nmcmc <- 1000000
nchains <- 100


history <- foreach(ichain = 1:nchains, .combine = cbind) %dorng% {
  history_onechain <- rep(0, nmcmc)
  state <- rinit()
  for (imcmc in 1:nmcmc){
    state <- MALA(state)
    history_onechain[imcmc] <- state$position[506]
  }
  history_onechain
}

history.df <- reshape2::melt(history)
names(history.df) <- c("iteration", "chain", "value")
history.df$chain <- rep(1:nchains, each = nmcmc)

save(nmcmc, nchains, stepsize, Lmax, history.df, file = file.path(filepath, "randomeffects.mala.longrun.RData"))
load(file = file.path(filepath, "randomeffects.mala.longrun.RData"))

tail(history.df)

library(tidyverse)
ggplot(history.df %>% filter(chain <= 10), aes(x = iteration, y = value, color = chain)) + geom_line() + theme_minimal() +
  theme(legend.position = "none")
# 