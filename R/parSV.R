
##FFT function for calculating RSVe on a centered matrix A
mSVEfft <- function (A, b, method = "bartlett")
{
  library(fftwtools)
  n <- nrow(A) # A must be centered matrix
  p <- ncol(A)
  w <- as.vector(lag2(1:b, n = n, b = b, method = method)) # calculate lags
  w <- c(1, w[1:(n-1)], 0, w[(n-1):1])  # starting steps from FFT paper
  w <- Re(fftw_r2c(w))
  FF <- matrix(0, ncol = p, nrow = 2*n)
  FF[1:n,] <- A
  if(p > 1)  # multivariate
  {
    FF <- mvfftw_r2c (FF)
    FF <- FF * matrix(w, nrow = 2*n, ncol = p)
    FF <- mvfftw_c2r(FF) / (2* n )
    return ((t(A) %*% FF[1:n, ]) / n )
  } else if(p == 1)  ##univariate calls
  {
    FF <- fftw_r2c (FF)
    FF <- FF * matrix(w, nrow = 2*n, ncol = p)
    FF <- fftw_c2r(FF) / (2* n )
    return ((t(A) %*% FF[1:n]) / n )
  }
  
}

#'@export
parSVE <- function(chains, r = 1)
{
  # number of chains
  nchains <- length(chains)
  
  # finding batch size
  b.final <- 0
  for(m in 1:nchains)
  {
    b.final <- b.final + batchSize(chains[[m]])
  }
  b.final <- floor(b.final/nchains)
  
  global.mean <- mean(sapply(chains, mean))
  n <- length(chains)
  
  rsve <- 0
  
  for (m in 1:nchains)
  {
    chain.cen <- scale(chains[[m]], center = global.mean, scale =FALSE)
    foo <- mSVEfft(A = chain.cen, b = b.final, method = "tukey") # change to  "bartlett"
    
    rsve <- rsve + (2*foo - mSVEfft(A = chain.cen, b = floor(b.final/r), method = "tukey"))
  }
  
  rtn <- rsve/nchains
  return(rtn)
}
