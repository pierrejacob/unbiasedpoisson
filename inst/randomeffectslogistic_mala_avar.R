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



k <- 10000
lag <- k
m <- 5*k

# res <- sample_unbiasedvar_reservoir(MALA, CoupledMALA, rinit,  h = h, k = k, m = m, lag = lag, x_0 = x_0, natoms = 10)
# # names(res)
# # res$elapsedtime
# res$cost
# res$cost_umcmc
# res$cost_fishyterms
# # res$estimator

nrep <- 1e3
natoms_seq <- c(1, 5, 10)
results_ <-  foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  run <- sample_unbiasedvar_reservoir(MALA, CoupledMALA, rinit,  h = h, k = k, m = m, lag = lag, x_0 = x_0, natoms = natoms_seq)
  data.frame(sampler = "MALA",
             natoms = natoms_seq,
             rep = irep,
             estimator = run$estimator[1,],
             cost = run$cost,
             pih = mean(run$pih),
             varh = run$varh,
             cost_fishyestimation = run$cost_fishyterms)
}

save(results_, nrep, natoms_seq, k, m, lag, file = "output/randomeffects.mala.avar.RData")
load(file = "output/randomeffects.mala.avar.RData")
tail(results_)

## compute inefficiency of UMCMC
results_umcmc <- results_ %>% filter(natoms == 1) %>% mutate(cost_umcmc = cost - cost_fishyestimation) %>% select(pih, cost_umcmc)
cat("inefficiency of UMCMC:", round(var(results_umcmc$pih) * mean(results_umcmc$cost_umcmc), 3), "\n")
avar_bestguess <- results_ %>% filter(natoms == 10) %>% summarise(avar = mean(estimator)) %>% pull()
cat("inefficiency of std MCMC:", round(avar_bestguess, 3), "\n")

results_ %>% group_by(natoms) %>% summarise(meancost = mean(cost), varestimator = var(estimator)) %>% mutate(inef = meancost * varestimator)

###
results_ %>% summarise(meanpih = mean(pih), stderror = 2*sd(pih)/sqrt(nrep))

###
results_ %>% summarise(meanvarh = mean(varh), stderror = 2*sd(varh)/sqrt(nrep))

###
results_ %>% group_by(natoms) %>% summarise(meanavar = mean(estimator), varavar = var(estimator), stderror = 2*sd(estimator)/sqrt(nrep), meancost = mean(cost)) %>%
  mutate(inefficiency = format(varavar * meancost, digits = 2)) %>% select(-varavar)



get_mean_with_booterror <- function(v, nb = 1000){
  bootresult <- boot::boot(data = v, statistic = function(v,indices) mean(v[indices]), R = nb)
  paste0("[", 
         paste0(prettyNum(as.numeric(quantile(bootresult$t, probs = c(0.025, 0.975))), digits=2, scientific = F, nsmall = 0), collapse = " - "),
         "]")
}
# print variance
get_var_with_booterror <- function(v, nb = 1000){
  bootresult <- boot::boot(data = v, statistic = function(v,indices) var(v[indices]), R = nb)
  paste0("[", 
         paste0(prettyNum(as.numeric(quantile(bootresult$t, probs = c(0.025, 0.975))), digits=2, scientific = T), collapse = " - "),
         "]")
}
# print efficiency
get_eff_with_booterror <- function(v, c, nb = 1000){
  bootresult <- boot::boot(data = cbind(v,c), statistic = function(x,indices) var(x[indices,1])*mean(x[indices,2]), R = nb)
  paste0("[", 
         paste0(prettyNum(as.numeric(quantile(bootresult$t, probs = c(0.025, 0.975))), digits=2, scientific = T), collapse = " - "),
         "]")
}


library(dplyr)
table_ <- results_ %>% ungroup() %>% group_by(natoms) %>%
  summarise(pestimate = get_mean_with_booterror(estimator), pcost = get_mean_with_booterror(cost), pfishycost = get_mean_with_booterror(cost_fishyestimation), pvariance = get_var_with_booterror(estimator), pefficiency = get_eff_with_booterror(estimator,cost))
table <- table_ %>% setNames(c("R", "estimate", "total cost", "fishy cost", 
                               "variance of estimator", "inefficiency"))
print(table_, width = Inf)
# knitr::kable(table_, digits = 0, row.names = NA, format = 'latex', escape = FALSE) %>%
#   cat(., file = 'output/randomeffects.table.tex')



