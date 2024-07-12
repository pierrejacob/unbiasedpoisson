rm(list = ls())
library(unbiasedpoisson)
library(doParallel)
library(doRNG)
library(tidyverse)
library(mcmcse) 

set.seed(1)
setmytheme()
registerDoParallel(cores = detectCores() - 2)

load(file = "~/Dropbox/UnbiasedPoissonNumerics/binomialssm1d.longrun.RData")
history.df <- history %>% rename(value = x1) %>% select(-targetpdf)
history.df %>% tail
## discard first iterations
burnin <- 1000
## group chains by 4
chainindices <- matrix(1:nchains, ncol = 4)
nrep <- nrow(chainindices)
pmmhclassicalvardf <- data.frame()
for (nmcmc_ in c(1e4, 4e4)){
  
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
  pmmhclassicalvardf <- rbind(pmmhclassicalvardf, rbind(bmdf, svdf))
}
save(nmcmc, nchains, burnin, pmmhclassicalvardf, 
     file = "output/binomialssm.classicalavar.RData")
load(file = "output/binomialssm.classicalavar.RData")

pmmhclassicalvardf %>% group_by(nmcmc, method, r) %>% summarise(mean = mean(value), 
                                                                sd = sd(value),
                                                                v = var(value))

load("output/binomialssm1d.uavar.P256.goodanchor.RData")

res_df_goodanchor <- foreach (irep = 1:nrep, .combine = rbind) %do% {
  res_ <- results.goodx0[[irep]]
  run <- unbiasedvar_from_tworuns(res_[[1]], res_[[2]], 
                                  natoms = natoms_seq)
  data.frame(natoms = natoms_seq,
             rep = irep,
             estimator = run$estimator[1,],
             cost = run$cost,
             pih = mean(run$pih),
             varh = run$varh,
             cost_fishyestimation = run$cost_fishyterms)
}

res_df_goodanchor %>% tail
results.df <- res_df_goodanchor
results_uavar <- results.df %>% filter(natoms == 50) %>% pull(estimator)
estim_lo <- mean(results_uavar) - 1.96 * sd(results_uavar) / sqrt(length(results_uavar))
estim_hi <- mean(results_uavar) + 1.96 * sd(results_uavar) / sqrt(length(results_uavar))

g <- ggplot(pmmhclassicalvardf %>% filter(nmcmc == 10000) %>% mutate(tss = paste(method, r)),
            aes(x = tss, y = value)) + geom_point() +
  stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 3,
    shape = 24,
    fill = "yellow"
  )
g <- g + xlab("method") + ylab("estimate") + ylim(0.001, 0.004)
g <- g + geom_hline(yintercept =  mean(results_uavar), col = "red") +  geom_hline(yintercept = c(estim_lo, estim_hi), linetype = "dashed", color = "red")
g

g <- ggplot(pmmhclassicalvardf %>% filter(nmcmc == 40000) %>% mutate(tss = paste(method, r)),
            aes(x = tss, y = value)) + geom_point() +
  stat_summary(
    geom = "point",
    fun.y = "mean",
    col = "black",
    size = 3,
    shape = 24,
    fill = "yellow"
  )
g <- g + xlab("method") + ylab("estimate")  + ylim(0.001, 0.004)
g <- g + geom_hline(yintercept =  mean(results_uavar), col = "red") +  geom_hline(yintercept = c(estim_lo, estim_hi), linetype = "dashed", color = "red")
g

