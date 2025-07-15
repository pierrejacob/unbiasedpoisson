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
# source("inst/mala_functions.R")
# filepath <- "~/Dropbox/UnbiasedPoissonNumerics"
filepath <- ""

load(file = file.path(filepath, "randomeffects.mala.longrun.RData"))
tail(history.df)

## discard first iterations
burnin <- 15000
## group chains by 4
chainindices <- matrix(1:nchains, ncol = 4)
nrep <- nrow(chainindices)
malaclassicalvardf <- data.frame()
for (nmcmc_ in c(1e5, 1e6)){
  
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
  malaclassicalvardf <- rbind(malaclassicalvardf, rbind(bmdf, svdf))
}
save(nmcmc, nchains, burnin, malaclassicalvardf, 
     file = "output/randomeffects.mala.classicalavar.RData")
load(file = "output/randomeffects.mala.classicalavar.RData")

malaclassicalvardf %>% group_by(nmcmc, method, r) %>% summarise(mean = mean(value), 
                                                               sd = sd(value),
                                                               v = var(value))


load(file = "output/randomeffects.mala.avar.RData")

results.mala <- results_ %>% filter(natoms == 10) %>% pull(estimator)
estim_lo <- mean(results.mala) - 1.96 * sd(results.mala) / sqrt(length(results.mala))
estim_hi <- mean(results.mala) + 1.96 * sd(results.mala) / sqrt(length(results.mala))

g <- ggplot(malaclassicalvardf %>% filter(nmcmc == 100000) %>% mutate(tss = paste(method, r)),
            aes(x = tss, y = value)) + geom_point() +
  stat_summary(
    geom = "point",
    fun.y = "mean",
    col = "black",
    size = 3,
    shape = 24,
    fill = "yellow"
  )
g <- g + xlab("method") + ylab("estimate") + ylim(50, 150)
g <- g + geom_hline(yintercept =  mean(results.mala), col = "red") +  geom_hline(yintercept = c(estim_lo, estim_hi), linetype = "dashed", color = "red")
g

g <- ggplot(malaclassicalvardf %>% filter(nmcmc == 1000000) %>% mutate(tss = paste(method, r)),
            aes(x = tss, y = value)) + geom_point() +
  stat_summary(
    geom = "point",
    fun.y = "mean",
    col = "black",
    size = 3,
    shape = 24,
    fill = "yellow"
  )
g <- g + xlab("method") + ylab("estimate") + ylim(50, 150)
g <- g + geom_hline(yintercept =  mean(results.mala), col = "red") +  geom_hline(yintercept = c(estim_lo, estim_hi), linetype = "dashed", color = "red")
g

