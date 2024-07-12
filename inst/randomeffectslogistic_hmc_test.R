rm(list = ls())
library(unbiasedpoisson)
library(doParallel)
library(doRNG)
set.seed(1)
setmytheme()
registerDoParallel(cores = detectCores() - 2)
source("inst/randomeffectslogistic_model.R")
source("inst/hmc_functions.R")


### run HMC
# state <- rinit()
# nmcmc <- 1000
# betahistory <- rep(0, nmcmc)
# betahistory[1] <- state$position[506]
# naccepts <- 0
# for (i in 2:nmcmc){
#   state <- HMC(state)
#   naccepts <- naccepts + state$accept
#   betahistory[i] <- h(state$position)
#   if (i %% 100 == 0){
#     cat("Iteration ", i, "\n")
#     cat("Accept rate:", naccepts/i*100, "%\n")
#   }
# }
# matplot(betahistory, type = "l", lty = 1, col = 1, xlab = "Iteration", ylab = "beta[506]")


state1 <- rinit()
state2 <- rinit()
nmcmc <- 300
betahistory <- matrix(0, nrow = nmcmc, ncol = 2)
for (imcmc in 1:nmcmc){
  res_ <- CoupledHMC(state1, state2)
  state1 <- res_$state1
  state2 <- res_$state2
  betahistory[imcmc,1] <- state1$position[506]
  betahistory[imcmc,2] <- state2$position[506]
}
matplot(betahistory[,1] - betahistory[,2], type = 'l')

sample_meetingtime(HMC, CoupledHMC, rinit, max_iterations = 1000)

# # meeting times
# lag <- 300
# nrep <- 500
# mts <- foreach(irep = 1:nrep, .combine = c) %dorng% {
#   res <- sample_meetingtime(HMC, CoupledHMC, rinit, lag = lag)
#   res$meetingtime
# }
# hist(mts)
# mean(mts-lag)
# ### note: sometimes error:
# ### Error in if (accept2) { : missing value where TRUE/FALSE needed
# 
# ## plot TV upper bounds
# niterations <- max(mts)
# ubounds <- sapply(1:niterations, function(t) tv_upper_bound(mts, lag, t))
# g_tvbounds <- qplot(x = 1:niterations, y = ubounds, geom = "line") +
#   ylab("TV distance") + xlab("t") + geom_hline(yintercept = 1)
# g_tvbounds <- g_tvbounds + scale_y_log10(breaks = c(0, 0.001, 0.01, 0.1,1)) + geom_rangeframe()
# g_tvbounds


k <- 500
lag <- 500
m <- 5*k
state0 <- rinit()
x0 <- state0$position
h <- function(x) x[506]

# res <- sample_unbiasedvar_reservoir(HMC, CoupledHMC, rinit,  h = h, k = k, m = m, lag = lag, x_0 = x0, natoms = 10)
# names(res)
# res$elapsedtime
# res$cost
# res$cost_umcmc
# res$cost_fishyterms
# res$estimator

nrep <- 1e3
natoms_seq <- c(1, 5, 10)
results_ <-  foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  run <- sample_unbiasedvar_reservoir(HMC, CoupledHMC, rinit,  h = h, k = k, m = m, lag = lag, x_0 = x0, natoms = natoms_seq)
  data.frame(sampler = "HMC",
             natoms = natoms_seq,
             rep = irep,
             estimator = run$estimator[1,],
             cost = run$cost,
             pih = mean(run$pih),
             varh = run$varh,
             cost_fishyestimation = run$cost_fishyterms)
}

save(results_, nrep, natoms_seq, k, m, lag, file = "output/randomeffects.hmc.avar.RData")
load(file = "output/randomeffects.hmc.avar.RData")
head(results_)

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



