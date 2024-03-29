% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/highdim_eta_update.R
\name{crn_eta_update}
\alias{crn_eta_update}
\title{eta_update_half_t_crn_couple}
\usage{
crn_eta_update(
  xi_1,
  Beta_1,
  eta_1,
  sigma2_1,
  xi_2,
  Beta_2,
  eta_2,
  sigma2_2,
  t_dist_df
)
}
\arguments{
\item{xi_1, xi_2}{xi values from the pair of chains}

\item{Beta_1, Beta_2}{beta values (each vector of length p) from the pair of chains}

\item{eta_1, eta_2}{eta values (each vector of length p) from the pair of chains}

\item{sigma2_1, sigma2_2}{sigma values from the pair of chains}

\item{t_dist_df}{half-t degree of freedom}
}
\value{
Returns (eta_1, eta_2) under common random numbers coupling
}
\description{
Common random numbers coupling of Eta Update
}
