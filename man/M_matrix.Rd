% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/highdim_kernel.R
\name{M_matrix}
\alias{M_matrix}
\title{M_matrix}
\usage{
M_matrix(xi, eta, X_eta_tX_matrix, n)
}
\arguments{
\item{xi}{positive scalar}

\item{eta}{vector of length p}

\item{X_eta_tX_matrix}{X * Diag(1/eta) * t(X_transpose), matrix n by n}

\item{n}{positive integer}
}
\value{
I + X Diag(1/eta) t(X_transpose)
}
\description{
Calculates I + X * Diag(1/eta) * t(X_transpose)
}
