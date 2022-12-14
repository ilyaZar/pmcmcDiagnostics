% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02_table_diagnostics.R
\name{compute_ess}
\alias{compute_ess}
\title{Computes the expected sample size.}
\usage{
compute_ess(mcmc_sims_after, ESS_STANDARD, ESS_STAN)
}
\arguments{
\item{mcmc_sims_after}{(particle) MCMC draws/simulated samples after burn-in}

\item{ESS_STANDARD}{logical; if \code{TRUE} computes the effective sample
size}

\item{ESS_STAN}{logical; if \code{TRUE} computes the effective sample
size in Stan-style i.e. with ESS bulk and ESS tail}
}
\value{
a vector of expected sample sizes of length equal to the number of
parameters
}
\description{
Uses samples after burn-in period.
}
