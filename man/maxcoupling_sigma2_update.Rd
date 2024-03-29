% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/highdim_coupledkernel.R
\name{maxcoupling_sigma2_update}
\alias{maxcoupling_sigma2_update}
\title{xi_update}
\usage{
maxcoupling_sigma2_update(xi_1, eta_1, xi_2, eta_2, n, ssr_1, ssr_2, a0, b0)
}
\arguments{
\item{xi_1, xi_2}{current xi values (positive scalars)}

\item{eta_1, eta_2}{current eta values (vector length p)}

\item{n}{number of data points}

\item{ssr_1, ssr_2}{postiive scalars}

\item{a0}{positive scalar}

\item{b0}{positive scalar}
}
\value{
Coupled update of sigma2 given eta
}
\description{
Coupled update of sigma2 given eta
}
