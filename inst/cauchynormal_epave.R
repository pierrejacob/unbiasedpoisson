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

## import functions for the Cauchy-Normal model
source("inst/cauchynormal_functions.R")
## 
load(file = "output/cauchynormal.mrth.meetings.RData")
niterations <- max(meetingtimes_mrth-lag)
ubounds_mrth <- sapply(1:niterations, function(t) tv_upper_bound(unlist(meetingtimes_mrth), lag, t))
g_tvbounds <- qplot(x = 1:niterations, y = ubounds_mrth, geom = "line") +
  ylab("TV distance") + xlab("t")
g_tvbounds <- g_tvbounds + scale_y_log10(breaks = c(0, 0.001, 0.01, 0.1,1)) + geom_rangeframe() + geom_hline(yintercept = 1)
g_tvbounds

## 
cat("mean meeting time:", round(mean(meetingtimes_mrth-lag), 2), "\n")

# state <- rinit_mrth()
# state <- single_kernel_mrth(state)
# coupled_kernel_mrth(state, rinit_mrth())



nmcmc <- 10000
burnin <- 100
## tuning parameter D:
## estimate fishy function only every D iterations
## cost : nmcmc for the MCMC part
## and 2 * (nmcmc/D) * E[tau] for the G part
## we want to balance the two costs: nmcmc = 2 * (nmcmc/D) * E[tau]
## i.e. D = 2 * E[tau]
D <- ceiling(mean(meetingtimes_mrth-lag)) * 2

nrep <- 1000
results_epave <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  run <- sample_epave(single_kernel_mrth, coupled_kernel_mrth, rinit_mrth, 
               nmcmc, burnin, D, h, state_x_0_mrth)
  data.frame(sampler = "mrth",
             D = D,
             rep = irep,
             estimator = run$vEPAVE[1,1],
             cost = run$costEPAVE,
             pih = run$meanh,
             varh = run$vMC,
             cost_fishyestimation = run$costG)
}
save(results_epave, nrep, D, nmcmc, burnin, file = "output/cauchynormal.mrth.epave.RData")
load(file = "output/cauchynormal.mrth.epave.RData")

results_epave %>% tail
results_epave %>% summarise(avar = mean(estimator), 
                            varavar = var(estimator), varh = mean(varh), meancost = mean(cost))


### verification with SUAVE
load(file = "output/cauchynormal.mrth.avar.RData")
results.mrth %>% tail
results.mrth %>% filter(natoms == 100) %>% summarise(avar = mean(estimator), 
                                                     varavar = var(estimator),
                                                     varh = mean(varh),
                                                     meancost = mean(cost))

