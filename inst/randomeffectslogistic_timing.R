rm(list = ls())
library(unbiasedpoisson)
library(doParallel)
library(doRNG)
set.seed(1)
setmytheme()
registerDoParallel(cores = detectCores() - 2)
source("inst/randomeffectslogistic_model.R")
source("inst/hmc_functions.R")

nmcmc <- 100
nrep <- 10
timingHMC <- c()
for (rep in 1:nrep){
  ## start timer
  start_time <- Sys.time()
  state <- rinit()
  for (imcmc in 1:nmcmc){
    state <- HMC(state)
  }
  ## record elapsed time in seconds
  timingHMC <- c(timingHMC, as.numeric(Sys.time() - start_time, units = "secs"))
}

summary(timingHMC/nmcmc)


nmcmc <- 100
nrep <- 10
timingMALA <- c()
for (rep in 1:nrep){
  ## start timer
  start_time <- Sys.time()
  state <- rinit()
  for (imcmc in 1:nmcmc){
    state <- MALA(state)
  }
  ## record elapsed time in seconds
  timingMALA <- c(timingMALA, as.numeric(Sys.time() - start_time, units = "secs"))
}

summary(timingMALA/nmcmc)

mean(timingMALA)
mean(timingHMC)

## time gradtarget
nrep <- 100
timing_gradtarget <- c()
for (rep in 1:nrep){
  ## start timer
  state <- rinit()
  start_time <- Sys.time()
  ## 
  gradtarget(state$position)
  ## record elapsed time in seconds
  timing_gradtarget <- c(timing_gradtarget, as.numeric(Sys.time() - start_time, units = "secs"))
}

summary(timing_gradtarget)



