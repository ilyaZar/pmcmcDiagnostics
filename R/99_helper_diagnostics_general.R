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
#' @inheritParams get_update_rates
#'
#' @return a matrix of update rates of appropriate dimensions: TT x NN
#' @export
compute_states_urs <- function(
  trajectories,
  dim_list = list("TT" = 1, "DD" = 2, "MM" = 3, "NN" = 4),
  WITH_CHECKS = FALSE) {
  dims_tkn <- dim(trajectories)
  TT <- dims_tkn[dim_list[["TT"]]]
  MM <- dims_tkn[dim_list[["MM"]]]
  NN <- dims_tkn[dim_list[["NN"]]]
  trajectories <- aperm(trajectories, c(dim_list["TT"],
                                        dim_list["DD"],
                                        dim_list["MM"],
                                        dim_list["NN"]))
  urs_all <- urs_big(trajectories, TT, NN, MM, WITH_CHECKS)
  return(urs_all)
}
adjust_trajectories <- function(traj_mat) {
  ID_ZERO_CMPNTS_MAT <- apply(traj_mat, c(1, 2), check_zero_component_traj)
  stopifnot(`FAILED.` = check_zero_cmpnt_mat_fail(ID_ZERO_CMPNTS_MAT))
  id_tkn <- 1
  while (id_tkn <= nrow(ID_ZERO_CMPNTS_MAT)) {
    ID_ZERO_CMPNTS <- ID_ZERO_CMPNTS_MAT[id_tkn, ]
    if (all(ID_ZERO_CMPNTS)) {
      id_tkn <- id_tkn + 1
    } else {
      break
    }
  }
  if (id_tkn == nrow(ID_ZERO_CMPNTS_MAT)) stop("No variation in states ...")
  traj_mat[, !ID_ZERO_CMPNTS, ]
}
check_zero_cmpnt_mat_fail <- function(x) {
  rws       <- rowSums(x)
  rws       <- setdiff(rws, ncol(x))
  rws_unq   <- unique(rws)
  CHECK_LEN <- length(rws_unq) == 1
  if (CHECK_LEN) return(TRUE)
  FALSE
}
check_zero_component_traj <- function(x) {
  if (length(unique(x)) == 1) return(TRUE)
  return(FALSE)
}
urs_big <- function(trajectories, TT, NN, MM, WITH_CHECK = FALSE) {
  urs_all   <- matrix(0, nrow = TT, ncol = NN)
  for (n in 1:NN) {
    tmp_traj  <- adjust_trajectories(trajectories[, , , n])
    if (WITH_CHECK) {
      DD <- ncol(tmp_traj)
      urs_all[, n] <- compute_urs_with_checks(tmp_traj, TT, DD, MM)
    } else {
      urs_all[, n] <- urs_small(tmp_traj[, 1, ], MM)
    }
  }
  return(urs_all)
}
urs_small <- function(traj_mat, num_draws) {
  apply(traj_mat, MARGIN = 1, compute_frac_unq, num_draws, simplify = TRUE)
}
compute_frac_unq <- function(x, num_draws) {
  length(unique(x)) / num_draws
}
compute_urs_with_checks <- function(x, TT, DD, MM) {
  urs_per_n <- matrix(0, nrow = TT, ncol = DD)
  for (dd in seq_len(DD)) {
    urs_per_n[, dd] <- urs_small(x[, dd, ], MM)
  }
  unique_states <- apply(urs_per_n,
                         MARGIN = 1,
                         function(x) {
                           Reduce(equal_values, x)
                         })
  CHECK_COMPUTATIONS <- any(unique_states != urs_per_n[, 1])
  if (CHECK_COMPUTATIONS) {
    browser()
    stop("Numerical problems during update rate computations!")
  }
  return(unique_states)
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
  if (is.null(list_true_vals) || all(is.na(list_true_vals))) return(NA)
  names_pars <- names(list_true_vals)
  CHECK_DIST_SPECIAL <- grepl("Gen", class(list_true_vals)[1])
  if (CHECK_DIST_SPECIAL) {
    out <- get_true_vals_special(list_true_vals, names_pars)
  } else {
    out <- get_true_vals_default(list_true_vals, names_pars)
  }
  return(out)
}
get_true_vals_special <- function(list_true_vals, names_pars) {
  DD <- nrow(list_true_vals[["sig_sq"]])
  out <- numeric(0)
  for(j in names_pars) {
    for (d in 1:DD) {
      if (is.null(list_true_vals[[j]])) next;
      if (j == c("sig_sq")) {
        out <- c(out, list_true_vals[[j]][d, 1, 1])
        out <- c(out, list_true_vals[[j]][d, 1, 2])
      } else if (j == "phi") {
        out <- c(out, unlist(list_true_vals[[j]][["A"]][[d]][, 1]))
        out <- c(out, unlist(list_true_vals[[j]][["B"]][[d]][, 1]))
      } else if (j == "beta_z_lin") {
        out <- c(out, unlist(list_true_vals[[j]][["A"]][[d]]))
        out <- c(out, unlist(list_true_vals[[j]][["B"]][[d]]))
      } else if (j == "beta_u_lin") {
        out <- c(out, unlist(t(list_true_vals[[j]][["A"]][[d]])))
        out <- c(out, unlist(t(list_true_vals[[j]][["B"]][[d]])))
      } else if (j == "vcm_u_lin") {
        out <- c(out, unlist(list_true_vals[[j]][["A"]][[d]]))
        out <- c(out, unlist(list_true_vals[[j]][["B"]][[d]]))
      }
    }
  }
  out <- unname(out)
  return(out)
}
get_true_vals_default <- function(list_true_vals, names_pars) {
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
