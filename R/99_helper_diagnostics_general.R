#' Computation of update rates for particle states
#'
#' Computes the update rates for simulated/drawn particle trajectories of the
#' particle MCMC output. Update rates are a measure of how frequently the
#' particular state (as displayed along the x axes) gets updated in the MCMC
#' draws of the particle Gibbs/MH approach. The more frequently, the better as
#' higher update rates indicate a better mixing of the sampler.
#'
#' @param trajectories latent state trajectories (with ordering of mcmc draws
#'   and time dimension specified under states_in_cols)
#' @param states_in_cols logical with default TRUE; if TRUE, mcmc draws are
#'   given in rows and time dimension (i.e. number of states) in columns
#' @param return_values logical; if TRUE, the update rates are returned as
#'   values rather than a plot
#'
#' @return object of type matplot displaying the udpate rates of particle
#'   trajectories
#' @export
analyse_states_ur <- function(trajectories,
                              states_in_cols = TRUE,
                              return_values = FALSE) {
  dim_trajs <- length(dim(trajectories))
  if (dim_trajs == 4) {
    NN <- dim(trajectories)[4]
  } else {
    NN <- 1
  }
  if (states_in_cols) {
    num_states <- dim(trajectories)[1]
    num_comps <-  dim(trajectories)[2]
    num_draws <-  dim(trajectories)[3]
  } else if (!states_in_cols) {
    num_states <- dim(trajectories)[2]
    num_comps <-  dim(trajectories)[1]
    num_draws <-  dim(trajectories)[3]
    if (dim_trajs == 4) {
      trajectories <- aperm(trajectories, c(2, 1, 3, 4))
    } else {
      trajectories <- aperm(trajectories, c(2, 1, 3))
    }
  }

  urs_all   <- matrix(0, nrow = num_states, ncol = NN)

  for (n in 1:NN) {
    urs_per_n <- matrix(0, nrow = num_states, ncol = num_comps)
    for (d in 1:num_comps) {
      num_unique_states <- apply(trajectories[ ,d , ,n], MARGIN = 1, unique)
      #apply(trajectories[ ,d , ,n], MARGIN = 2, unique)
      num_unique_states <- unlist(lapply(num_unique_states, length))
      urs_per_n[, d] <- num_unique_states/num_draws
    }
    unique_states <- apply(urs_per_n, MARGIN = 1,
                               function(x) {Reduce(equal_values, x)})
    test_computations <- any(unique_states != urs_per_n[, 1])
    # browser()
    if (test_computations) {
      stop("Numerical problems during update rate computations!")
    } else {
      urs_all[, n] <- unique_states
    }
  }

  # browser()
  # urs[1, ] <- min(urs[2:num_states, 1])
  # .colMeans(m = num_states, n = num_comps)
  if (return_values) {
    return(urs_all)
  }
  graphics::matplot(urs_all , type = "l")
}
#' To be used within \code{Reduce()}-constructs to check if values are equal
#'
#' Equal values as returned by \code{all.equal()}, not \code{identical()}.
#'
#' @param x first value
#' @param y second value (that is checked to be equal to \code{x})
#'
#' @return either \code{TRUE} if \code{x} and \code{y} are equal or \code{FALSE}
#'   otherwise
equal_values <- function(x ,y) {
  if (isTRUE(all.equal(x, y, check.attributes = FALSE))) {
    return(x)
  } else {
    return(FALSE)
  }
}
#' Eliminate burnin and apply thinning to parameter draws of MCMC output
#'
#' Given a matrix of mcmc draws, returns a matrix with draws after burnin period
#' and thinning.
#'
#' @param draws parameter draws of the (particle) MCMC approach in matrix form
#'   where the number of draws are stored in rows and the different parameters
#'   in columns
#' @param burnin burnin period (e.g. 1000 = the first 1000 mcmc draws are
#'   eliminated from the chain)
#' @param thin thinning frequency (e.g. 10 = every 10th mcmc draw is eliminated
#'   from the chain)
#'
#' @return a matrix of (particle) MCMC parameter draws after burnin and thinning
#' @export
burn_and_thin <- function(draws, burnin, thin = NULL) {
  # burnin <- burnin + 1
  burned_interval  <- burnin:nrow(draws)
  mcmc_sims_after  <- draws[burned_interval, ]
  if (!(is.null(thin))) {
    thinned_interval <- seq(from = 1, to = nrow(mcmc_sims_after), by = thin)
    mcmc_sims_after <- mcmc_sims_after[thinned_interval, ]
  }
  return(mcmc_sims_after)
}
