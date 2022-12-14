% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/99_helper_pgas_output_transformation.R
\name{pgas_out_2_list}
\alias{pgas_out_2_list}
\title{Transforms PGAS output to list format for comparison/testing reasons}
\usage{
pgas_out_2_list(pgas_out, DD, NN, TT, MM, dim_bet_z, cpp = NULL)
}
\arguments{
\item{pgas_out}{output of \code{pgas_cpp()} or \code{pgas_R()}}

\item{DD}{multivariate response/measurement dimension: e.g. number of
shares/fractions if the measurements are from a Dirichlet}

\item{NN}{cross sectional dimension}

\item{TT}{number of time periods}

\item{MM}{number of overall (particle) MCMC iterations}

\item{dim_bet_z}{number of regressors coefficients of z-type regressors for
each \code{d,...,DD}}

\item{cpp}{logical; if \code{TRUE}, the particle gibbs output comes from cpp
code, if \code{FALSE} from R code}
}
\value{
a list of 19 elements corresponding to former output format
}
\description{
Transforms PGAS output to list format for comparison/testing reasons
}
