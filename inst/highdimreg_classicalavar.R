rm(list = ls())
library(unbiasedpoisson)
library(doParallel)
library(doRNG)
library(tidyverse)
library(mcmcse) 
library(coda)

set.seed(1)
setmytheme()
registerDoParallel(cores = detectCores() - 2)

# filepath <- "~/Dropbox/UnbiasedPoissonNumerics"
filepath <- ""

load(file = file.path(filepath,"highdimreg.longrun.RData"))
history.df %>% head
## discard first iterations
burnin <- 1500


spectrum0_df <- data.frame()
for (nmcmc_ in c(2e4, 2e5)){
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
gibbsclassicalvardf <- data.frame()
for (nmcmc_ in c(2e4, 2e5)){
  
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
  gibbsclassicalvardf <- rbind(gibbsclassicalvardf, rbind(bmdf, svdf))
}


gibbsclassicalvardf <- gibbsclassicalvardf %>% mutate(tss = paste(method, r)) %>% select(value, nmcmc, tss)
gibbsclassicalvardf <- rbind(spectrum0_df %>% rename(tss = method) %>%  select(value, nmcmc, tss), gibbsclassicalvardf)
gibbsclassicalvardf$tss <- factor(gibbsclassicalvardf$tss, levels = c("spectrum0", "BM r = 1", "BM r = 2", "BM r = 3", "SV r = 1", "SV r = 2", "SV r = 3"))


save(nmcmc, nchains, burnin, gibbsclassicalvardf, 
     file = "output/highdimreg.classicalavar.RData")
load(file = "output/highdimreg.classicalavar.RData")

gibbsclassicalvardf %>% group_by(nmcmc, tss) %>% summarise(mean = mean(value), 
                                                                sd = sd(value),
                                                                v = var(value))

load(file = "output/highdimreg.uavar.RData")

results.df <- foreach(irep = 1:nrep, .combine=rbind) %do% {
  run <- results[[irep]]
  data.frame(natoms = natoms_seq,
             rep = irep,
             estimator = run$estimator[1,],
             cost = run$cost,
             pih = mean(run$pih),
             varh = run$varh,
             cost_fishyestimation = run$cost_fishyterms)
}

results_uavar <- results.df %>% filter(natoms == 10) %>% pull(estimator)
estim_lo <- mean(results_uavar) - 1.96 * sd(results_uavar) / sqrt(length(results_uavar))
estim_hi <- mean(results_uavar) + 1.96 * sd(results_uavar) / sqrt(length(results_uavar))

g <- ggplot(gibbsclassicalvardf %>% filter(nmcmc == 20000),
            aes(x = tss, y = value)) + geom_point() +
  stat_summary(
    geom = "point",
    fun.y = "mean",
    col = "black",
    size = 3,
    shape = 24,
    fill = "yellow"
  )
g <- g + xlab("method") + ylab("estimate") + ylim(25, 155)
g <- g + geom_hline(yintercept =  mean(results_uavar), col = "red") +  geom_hline(yintercept = c(estim_lo, estim_hi), linetype = "dashed", color = "red")
g

g <- ggplot(gibbsclassicalvardf %>% filter(nmcmc == 200000),
            aes(x = tss, y = value)) + geom_point() +
  stat_summary(
    geom = "point",
    fun.y = "mean",
    col = "black",
    size = 3,
    shape = 24,
    fill = "yellow"
  )
g <- g + xlab("method") + ylab("estimate") + ylim(25, 155)
g <- g + geom_hline(yintercept =  mean(results_uavar), col = "red") +  geom_hline(yintercept = c(estim_lo, estim_hi), linetype = "dashed", color = "red")
g

