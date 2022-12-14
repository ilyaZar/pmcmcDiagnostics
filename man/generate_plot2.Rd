% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/99_helper_diagnostics_mcmc_plot_helper.R
\name{generate_plot2}
\alias{generate_plot2}
\title{Generates base graphic analysis plot of MCMC output.}
\usage{
generate_plot2(
  mcmc_sims,
  mcmc_sims_after,
  burn = 0,
  thin = 0,
  num_mcmc,
  par_names,
  par_names_plots,
  true_vals = NULL,
  posterior_means,
  plot_num
)
}
\arguments{
\item{mcmc_sims}{mcmc draws as a matrix}

\item{mcmc_sims_after}{mcmc draws as a data.frame after burn-in and
thinning}

\item{burn}{burn-in period}

\item{thin}{thinning; \code{thin = 10} means every 10th draw is excluded}

\item{num_mcmc}{total number of MCMC draws}

\item{par_names}{parameter names as used in the pgas output (col names for
\code{mcmcm_sims} e.g.)}

\item{par_names_plots}{parameter names for plot labeling}

\item{true_vals}{true values if a simulation is run; default is NULL and true
values will not be added for the histogram and trace plots}

\item{posterior_means}{posterior means (have to be pre-computed and passed
directly here)}

\item{plot_num}{number of parameter or plot}
}
\value{
base plot in 2x2 form displaying the mcmc-diagnostics
}
\description{
Generates base graphic analysis plot of MCMC output for a model parameter
including: 1. histogram of mcmc draws (after burn-in and thinning) 2. trace
plot of mcmc draws (after burn-in and thinning) 3. autocorrelation plot of
mcmc draws (after burn-in and thinning) 4. trace plot of mcmc draws (BEFORE
burn-in and thinning)
}
