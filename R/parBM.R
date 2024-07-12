#'@rdname parBM
#'@title Parallel batch means
#'@description Parallel batch means
#' chain = list of chains, each of equal length
#'\itemize{
#' \item r = 1 is regular BM
#' \item r = 2 is flat-top BM
#' \item r = 3 is lugsail BM
#' }
#'@export
parBM <- function(chains, r = 1){
  library(mcmcse)
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
