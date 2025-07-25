rm(list = ls())
library(unbiasedpoisson)
setmytheme()
library(doParallel)
library(doRNG)
library(dplyr)
# register parallel cores
registerDoParallel(cores = 10)
# set RNG seed
set.seed(1)

source("inst/highdimreg_functions.R")
# filepath <- "~/Dropbox/UnbiasedPoissonNumerics"
filepath <- ""

nmcmc <- 2e5
nchains <- 100

history <- foreach(ichain = 1:nchains, .combine = cbind) %dorng% {
  history_onechain <- rep(0, nmcmc)
  state <- rinit()
  for (imcmc in 1:nmcmc){
    state <- single_kernel(state)
    history_onechain[imcmc] <- h(state$position)
  }
  history_onechain
}

history.df <- reshape2::melt(history)
names(history.df) <- c("iteration", "chain", "value")
history.df$chain <- rep(1:nchains, each = nmcmc)


save(nmcmc, nchains, history.df, file = file.path(filepath,"highdimreg.longrun.RData"))
load(file = file.path(filepath,"highdimreg.longrun.RData"))

library(mcmcse) 
##### Parallel batch means
# chain = list of chains, each of equal length
# r = 1 is regular BM
# r = 2 is flat-top BM
# r = 3 is lugsail BM
parBM <- function(chains, r = 1){
  # number of chains
  nchains <- length(chains)
  # finding batch size
  b.final <- 0
  for (m in 1:nchains){
    b.final <- b.final + batchSize(chains[[m]])
  }
  b.final <- floor(b.final/nchains)
  n <- length(chains[[1]])
  a <- floor(n/b.final)
  ab <- b.final * a
  trash <- n-ab
  big.chain <- numeric(length = ab*nchains)
  if (ab != n){
    for (i in 1:m){
      big.chain[((i-1)*ab+1):(i*ab)] <- chains[[i]][-(1:trash)]
    }
  } else {
    big.chain <- Reduce("c", chains)
  }
  rtn <- mcse.multi(big.chain, r = r)$cov
  return(rtn)
}

head(history.df)
dim(history.df)
history.df[nmcmc+1,]
burnin <- 1000
chainlist <- lapply(1:nchains, function(ichain) history.df[(nmcmc*(ichain-1)+burnin):(nmcmc*(ichain-1)+nmcmc),3])

parBM(chainlist, r = 1)
parBM(chainlist, r = 2)
parBM(chainlist, r = 3)

## BM for each chain
bms_ <- matrix(NA, nrow = nchains, ncol = 3)
for (ichain in 1:nchains){
  chainlist_ <- list(history.df[(nmcmc*(ichain-1)+burnin):(nmcmc*(ichain-1)+nmcmc),3])
  for (r in 1:3){
    bms_[ichain,r] <- parBM(chainlist_, r = r)
  }
}
hist(bms_[,1], prob = T, xlim = c(0, 100))
hist(bms_[,2], prob = T, col = "red", add = T)
hist(bms_[,3], prob = T, col = "blue", add = T)


# load(file = "output/highdimreg.uavar.RData")
# 
# results.df <- foreach(irep = 1:nrep, .combine=rbind) %do% {
#   run <- results[[irep]]
#   data.frame(natoms = natoms_seq,
#              rep = irep,
#              estimator = run$estimator[1,],
#              cost = run$cost,
#              pih = mean(run$pih),
#              varh = run$varh,
#              cost_fishyestimation = run$cost_fishyterms)
# }
# # 
# # hist(results.df %>% filter(natoms == 5) %>% pull(estimator), col = 'yellow', prob = T)

## variance of the target

var(history.df %>% filter(iteration > 100) %>% pull(value))
mean(history.df %>% filter(iteration > 100) %>% pull(value))

## test coupled kernel's validity 
# nmcmc <- 5e4
# nchains <- 10
# 
# # history <- foreach(ichain = 1:nchains, .combine = cbind) %dorng% {
# #   history_onechain <- rep(0, nmcmc)
# #   state <- rinit()
# #   for (imcmc in 1:nmcmc){
# #     state1 <- rinit()
# #     coupled_results <- coupled_kernel(state1, state)
# #     state <- coupled_results$state2
# #     history_onechain[imcmc] <- h(state$position)
# #   }
# #   history_onechain
# # }
# # 
# # history.df <- reshape2::melt(history)
# # names(history.df) <- c("iteration", "chain", "value")
# # history.df$chain <- rep(1:nchains, each = nmcmc)
# # save(nmcmc, nchains, history.df, file = "output/highdimreg.longrun_ckernel.RData")
# 
# load(file = "output/highdimreg.longrun_ckernel.RData")
# 
# load(file = "output/highdimreg.longrun.RData")
# history.df1 <- history.df
# hist(history.df %>% filter(iteration > 100) %>% pull(value), prob = T, col = rgb(1,0,1,0.5), nclass = 200)
# # acf(history.df %>% filter(iteration > 100, chain == 1) %>% pull(value))
# load(file = "output/highdimreg.longrun_ckernel.RData")
# history.df2 <- history.df
# hist(history.df %>% filter(iteration > 100) %>% pull(value), add = T, prob = T, col = rgb(1,1,0,0.5), nclass = 200)
# # acf(history.df %>% filter(iteration > 100, chain == 1) %>% pull(value))
# 
# ##
# tail(history.df1)
# tail(history.df2)
# ##
# acf(history.df1 %>% filter(iteration > 100, chain == 1) %>% pull(value))$acf[,,1]
# acf(history.df2 %>% filter(iteration > 100, chain == 1) %>% pull(value))$acf[,,1]
# ##
# mean(history.df1 %>% filter(iteration > 100) %>% pull(value) < -0.5)
# mean(history.df2 %>% filter(iteration > 100) %>% pull(value) < -0.5)
# ##
# 


