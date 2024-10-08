% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mvnorm_couplings.R
\name{rmvnorm_reflectionmax_diag}
\alias{rmvnorm_reflectionmax_diag}
\title{Reflection-Maximal coupling of two multivariate Normal distributions}
\usage{
rmvnorm_reflectionmax_diag(mu1, mu2, Cholesky)
}
\arguments{
\item{mu1}{First mean}

\item{mu2}{Second mean}

\item{Cholesky}{Cholesky factor, e.g. obtained with \code{\link[base]{chol}}}

\item{Cholesky_inverse}{Inverse of Cholesky factor, e.g. obtained with \code{solve(chol(Sigma))}}
}
\value{
A list containing 'xy', a matrix with 2 columns (one for each draw),
and a boolean indicator 'identical' indicating whether the two draws
are identical.
}
\description{
Sample from reflection-maximal coupling of two multivariate Normal distributions,
specified through their means, with the same covariance matrix, specified
through its Cholesky factor. Here the assumption is that the covariance is diagonal,
and thus the Cholesky factor is a vector of the same length as the means

The idea is that a multivariate Normal is drawn around the first mean (mu1),
and then reflected with respect to a hyperplane orthogonal to the direction between mu1 and mu2.

For univariate Normal distribution, see \code{\link{rnorm_reflectionmax}}.
}
\examples{
p <- 3
mu1 <- rep(0, p)
mu2 <- rep(1, p)
Sigma <- rep(0.4^2, p)
Sigma_chol <- sqrt(Sigma)
rmvnorm_reflectionmax_diag(mu1, mu2, Sigma_chol)
}
