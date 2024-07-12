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

load(file = "output/randomeffectslogistic.hmc.meetings.RData")

hist(mts-lag)


nmcmc <- 50000
burnin <- 500
## tuning parameter D:
## estimate fishy function only every D iterations
## cost : nmcmc for the MCMC part
## and 2 * (nmcmc/D) * E[tau] for the G part
## we want to balance the two costs: nmcmc = 2 * (nmcmc/D) * E[tau]
## i.e. D = 2 * E[tau]
D <- ceiling(mean(mts-lag)) * 2

nrep <- 100
results_epave <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  run <- sample_epave(HMC, CoupledHMC, rinit, 
                      nmcmc, burnin, D, h, state_x_0)
  data.frame(sampler = "HMC",
             D = D,
             rep = irep,
             estimator = run$vEPAVE,
             cost = run$costEPAVE,
             pih = run$meanh,
             varh = run$vMC,
             cost_fishyestimation = run$costG)
}
save(results_epave, nrep, D, nmcmc, burnin, file = "output/randomeffects.hmc.epave.RData")
load(file = "output/randomeffects.hmc.epave.RData")

results_epave %>% tail
results_epave %>% summarise(avar = mean(estimator), 
                            varavar = var(estimator), varh = mean(varh), meancost = mean(cost))

# load(file = "output/randomeffects.hmc.avar.RData")
# tail(results_)
# ## compute inefficiency of UMCMC
# results_umcmc <- results_ %>% filter(natoms == 1) %>% mutate(cost_umcmc = cost - cost_fishyestimation) %>% select(pih, cost_umcmc)
# cat("inefficiency of UMCMC:", round(var(results_umcmc$pih) * mean(results_umcmc$cost_umcmc), 3), "\n")
# avar_bestguess <- results_ %>% filter(natoms == 10) %>% summarise(avar = mean(estimator)) %>% pull()
# cat("inefficiency of std MCMC:", round(avar_bestguess, 3), "\n")
