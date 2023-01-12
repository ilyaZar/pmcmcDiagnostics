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

#' @return a matrix of update rates of appropriate dimensions
#' @export
compute_states_ur <- function(trajectories,
                              states_in_cols = TRUE) {
  dim_trajs <- length(dim(trajectories))
  if (dim_trajs == 4) {
    NN <- dim(trajectories)[4]
  } else {
    NN <- 1
  }
  if (states_in_cols) {
    num_states <- dim(trajectories)[1]
    # num_comps <-  dim(trajectories)[2]
    num_draws <-  dim(trajectories)[3]
  } else if (!states_in_cols) {
    num_states <- dim(trajectories)[2]
    # num_comps <-  dim(trajectories)[1]
    num_draws <-  dim(trajectories)[3]
    if (dim_trajs == 4) {
      trajectories <- aperm(trajectories, c(2, 1, 3, 4))
    } else {
      trajectories <- aperm(trajectories, c(2, 1, 3))
    }
  }
  urs_all <- urs_big(trajectories, NN, num_states, num_draws)
  return(urs_all)
}
adjust_trajectories <- function(traj_mat, num_draws) {
  check_components <- matrix(FALSE, nrow = num_draws,
                             ncol = ncol(traj_mat))
  for (mm in 1:num_draws) {
    check_components[mm, ] <- apply(traj_mat[, , mm, drop = FALSE], 2,
                                    function(x){length(unique(x)) == 1})
  }
  keep_cmp <- which(!apply(check_components, 2, function(x){all(x == TRUE)}))
  if (length(keep_cmp) == 0) stop("No variation in states; impoosible ...")
  traj_mat[, keep_cmp, ]
}
urs_big <- function(trajectories, NN, num_states, num_draws) {
  urs_all   <- matrix(0, nrow = num_states, ncol = NN)
  for (n in 1:NN) {
    tmp_traj  <- adjust_trajectories(trajectories[ , , , n],
                                     num_draws)
    num_comps <- dim(tmp_traj)[2]
    urs_per_n <- matrix(0, nrow = num_states, ncol = num_comps)
    for (d in 1:num_comps) {
      urs_per_n[, d] <- urs_small(tmp_traj[ , d, ], num_draws)
    }
    unique_states <- apply(urs_per_n, MARGIN = 1,
                           function(x) {Reduce(equal_values, x)})
    test_computations <- any(unique_states != urs_per_n[, 1])
    if (test_computations) {
      stop("Numerical problems during update rate computations!")
    } else {
      urs_all[, n] <- unique_states
    }
  }
  return(urs_all)
}
urs_small <- function(traj_mat, num_draws) {
  num_unique_states <- apply(traj_mat,
                             MARGIN = 1,
                             function(x) {length(unique(x))/ num_draws},
                             simplify = TRUE)
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
  burned_interval  <- burnin:nrow(draws)
  mcmc_sims_after  <- draws[burned_interval, ]
  if (!(is.null(thin))) {
    thinned_interval <- seq(from = 1, to = nrow(mcmc_sims_after), by = thin)
    mcmc_sims_after <- mcmc_sims_after[thinned_interval, ]
  }
  return(mcmc_sims_after)
}
get_true_vals <- function(list_true_vals) {
  if (is.null(list_true_vals) || is.na(list_true_vals)) return(NA)
  names_pars <- names(list_true_vals)

  DD <- nrow(list_true_vals[["sig_sq"]])
  out <- numeric(0)
  for(j in names_pars) {
    for (d in 1:DD) {
      if (is.null(list_true_vals[[j]])) next;
      if (j == c("sig_sq")) {
        out <- c(out, unique(list_true_vals[[j]][d, ]))
      } else if (j == "phi") {
        out <- c(out, unlist(list_true_vals[[j]][[d]][, 1]))
      } else if (j == "beta_z_lin") {
        out <- c(out, unlist(list_true_vals[[j]][[d]]))
      } else if (j == "beta_u_lin") {
        out <- c(out, unlist(t(list_true_vals[[j]][[d]])))
      } else if (j == "vcm_u_lin") {
        out <- c(out, unlist(list_true_vals[[j]][[d]]))
      }
    }
  }
  out <- unname(out)
  return(out)
}
