rm(list = ls())
library(unbiasedpoisson)
library(doParallel)
library(doRNG)
library(tidyverse)
library(mcmcse) 

set.seed(1)
setmytheme()
registerDoParallel(cores = detectCores() - 2)
# source("inst/randomeffectslogistic_model.R")
# source("inst/hmc_functions.R")
# filepath <- "~/Dropbox/UnbiasedPoissonNumerics"
filepath <- ""

load(file = file.path(filepath, "randomeffects.hmc.longrun.RData"))

tail(history.df)

## discard first iterations
burnin <- 1000



spectrum0_df <- data.frame()
for (nmcmc_ in c(1e4, 1e5)){
  spectrum0_df_ <- foreach (i = 1:nchains, .combine = rbind) %dopar% {
    xxx <- history.df %>% filter(iteration>burnin, iteration <= nmcmc_, chain == i) %>% pull(value)
    data.frame(method = "spectrum0", value = spectrum0(xxx)$spec)
  }
  spectrum0_df_$nmcmc <- nmcmc_
  spectrum0_df <- rbind(spectrum0_df, spectrum0_df_)
}

## group chains by 4
chainindices <- matrix(1:nchains, ncol = 4)
nrep <- nrow(chainindices)
hmcclassicalvardf <- data.frame()
for (nmcmc_ in c(1e4, 1e5)){
  
  bmdf <- foreach (irep = 1:nrep, .combine = rbind) %dopar% {
    chainlist <- lapply(chainindices[irep,], function(i) history.df %>% filter(iteration>burnin, iteration <= nmcmc_, chain == i) %>% pull(value))
    sapply(1:3, function(r) parBM(chainlist, r = r)[1,1])
  }
  bmdf <- as.data.frame(bmdf)
  svdf <- foreach (irep = 1:nrep, .combine = rbind) %dopar% {
    chainlist <- lapply(chainindices[irep,], function(i) history.df %>% filter(iteration>burnin, iteration <= nmcmc_, chain == i) %>% pull(value))
    sapply(1:3, function(r) parSVE(chainlist, r = r)[1,1])
  }
  svdf <- as.data.frame(svdf)
  names(bmdf) <- paste("r =", 1:3)
  names(svdf) <- paste("r =", 1:3)
  bmdf <- bmdf %>% gather(key = "r", value = "value") %>% mutate(method = "BM")
  bmdf$nmcmc <- nmcmc_
  svdf <- svdf %>% gather(key = "r", value = "value") %>% mutate(method = "SV")
  svdf$nmcmc <- nmcmc_
  hmcclassicalvardf <- rbind(hmcclassicalvardf, rbind(bmdf, svdf))
}

hmcclassicalvardf <- hmcclassicalvardf %>% mutate(tss = paste(method, r)) %>% select(value, nmcmc, tss)
hmcclassicalvardf <- rbind(spectrum0_df %>% rename(tss = method) %>%  select(value, nmcmc, tss), hmcclassicalvardf)
hmcclassicalvardf$tss <- factor(hmcclassicalvardf$tss, levels = c("spectrum0", "BM r = 1", "BM r = 2", "BM r = 3", "SV r = 1", "SV r = 2", "SV r = 3"))


save(nmcmc, nchains, burnin, hmcclassicalvardf,
     file = "output/randomeffects.hmc.classicalavar.RData")
load(file = "output/randomeffects.hmc.classicalavar.RData")

hmcclassicalvardf %>% group_by(nmcmc, tss) %>% summarise(mean = mean(value), 
                                                        sd = sd(value),
                                                        v = var(value))


load(file = "output/randomeffects.hmc.avar.RData")
results.hmc <- results_ %>% filter(natoms == 20) %>% pull(estimator)
estim_lo <- mean(results.hmc) - 1.96 * sd(results.hmc) / sqrt(length(results.hmc))
estim_hi <- mean(results.hmc) + 1.96 * sd(results.hmc) / sqrt(length(results.hmc))

g <- ggplot(hmcclassicalvardf %>% filter(nmcmc == 10000),
       aes(x = tss, y = value)) + geom_point() +
  stat_summary(
    geom = "point",
    fun.y = "mean",
    col = "black",
    size = 3,
    shape = 24,
    fill = "yellow"
  )
g <- g + xlab("method") + ylab("estimate") + ylim(1,2.8)
g <- g + geom_hline(yintercept =  mean(results.hmc), col = "red") +  geom_hline(yintercept = c(estim_lo, estim_hi), linetype = "dashed", color = "red")
g

g <- ggplot(hmcclassicalvardf %>% filter(nmcmc == 100000),
            aes(x = tss, y = value)) + geom_point() +
  stat_summary(
    geom = "point",
    fun.y = "mean",
    col = "black",
    size = 3,
    shape = 24,
    fill = "yellow"
  )
g <- g + xlab("method") + ylab("estimate")  + ylim(1,2.8)
g <- g + geom_hline(yintercept =  mean(results.hmc), col = "red") +  geom_hline(yintercept = c(estim_lo, estim_hi), linetype = "dashed", color = "red")
g
# table_hmc <- results.hmc %>% ungroup() %>% group_by(natoms) %>%
#   summarise(pestimate = get_mean_with_booterror(estimator), pcost = get_mean_with_booterror(cost), pfishycost = get_mean_with_booterror(cost_fishyestimation), pvariance = get_var_with_booterror(estimator), pefficiency = get_eff_with_booterror(estimator,cost))
# table_hmc <- table_hmc %>% setNames(c("R", "estimate", "total cost", "fishy cost", 
#                                       "variance of estimator", "inefficiency"))
# 


## what about spectrum0?
# spectrum0df <- data.frame()
# for (nmcmc_ in c(1e4, 1e5)){
#   spectrum0_ <- foreach (irep = 1:nchains, .combine = c) %dopar% {
#     coda::spectrum0(history.df %>% filter(iteration>burnin, iteration <= nmcmc_, chain == irep) %>% pull(value))$spec
#   }
#   spectrum0df <- rbind(spectrum0df, data.frame(nmcmc = nmcmc_, spectrum0 = spectrum0_, method = "spectrum0"))
# }
# names(spectrum0df) <- c("nmcmc", "value", "methodname")
# spectrum0df %>% head
# df_ <- rbind(hmcclassicalvardf %>% mutate(methodname = paste(method, r)) %>% select(nmcmc, value, methodname),
#              spectrum0df)
# g <- ggplot(df_,
#             aes(x = methodname, y = value)) + geom_point()
# g
## spectrum0 seems much worse than BM and SV