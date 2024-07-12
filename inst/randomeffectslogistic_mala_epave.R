rm(list = ls())
library(unbiasedpoisson)
library(doParallel)
library(doRNG)
library(tidyverse)
set.seed(1)
setmytheme()
registerDoParallel(cores = detectCores() - 2)
source("inst/randomeffectslogistic_model.R")
source("inst/mala_functions.R")

load(file = "output/randomeffectslogistic.mala.meetings.RData")

hist(mts-lag)


nmcmc <- 100000
burnin <- 15000
## tuning parameter D:
## estimate fishy function only every D iterations
## cost : nmcmc for the MCMC part
## and 2 * (nmcmc/D) * E[tau] for the G part
## we want to balance the two costs: nmcmc = 2 * (nmcmc/D) * E[tau]
## i.e. D = 2 * E[tau]
D <- ceiling(mean(mts-lag)) * 2

nrep <- 100
results_epave <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  run <- sample_epave(MALA, CoupledMALA, rinit, 
                      nmcmc, burnin, D, h, state_x_0)
  data.frame(sampler = "MALA",
             D = D,
             rep = irep,
             estimator = run$vEPAVE,
             cost = run$costEPAVE,
             pih = run$meanh,
             varh = run$vMC,
             cost_fishyestimation = run$costG)
}
save(results_epave, nrep, D, nmcmc, burnin, file = "output/randomeffects.mala.epave.RData")
load(file = "output/randomeffects.mala.epave.RData")

results_epave %>% tail
results_epave %>% summarise(avar = mean(estimator), 
                            varavar = var(estimator), varh = mean(varh), meancost = mean(cost))

# avar  varavar      varh meancost
# 115.3582 27977.93 0.1634143 233494.7

load(file = "output/randomeffects.mala.avar.RData")
tail(results_)
## 
results_ %>% filter(natoms == 10) %>% summarise(avar = mean(estimator), varavar = var(estimator), varh = mean(varh), meancost = mean(cost))
# avar  varavar      varh meancost
# 115.5854 21757.34 0.1645018 322786.9

## compute inefficiency of UMCMC
results_umcmc <- results_ %>% filter(natoms == 1) %>% mutate(cost_umcmc = cost - cost_fishyestimation) %>% select(pih, cost_umcmc)
cat("inefficiency of UMCMC:", round(var(results_umcmc$pih) * mean(results_umcmc$cost_umcmc), 3), "\n")
avar_bestguess <- results_ %>% filter(natoms == 10) %>% summarise(avar = mean(estimator)) %>% pull()
cat("inefficiency of std MCMC:", round(avar_bestguess, 3), "\n")

