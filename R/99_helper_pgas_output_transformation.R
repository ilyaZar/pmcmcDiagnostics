#' Transforms pgas output from \code{pgas_cpp()} or \code{pgas_R()} to a format
#' suitable for the function
#' \code{pmcmcDiagnostics::analyse_mcmc_convergence2()}
#'
#' @param pgas_out output from  \code{pgas_cpp()} or \code{pgas_R()}: a list of
#'   four elements:
#'   \itemize{
#'     \item{\code{sig_sq_xa:}}{matrix of dimension DD x num_mncmc_draws} of
#'     simulated standard deviation values
#'     \item{\code{phi_xa:}}{matrix of dimension DD x num_mncmc_draws} of
#'     simulated autoregressive parameter values
#'     \item{\code{bet_xa}}{matrix of dimension num_par_beta x num_mncmc_draws}
#'     of simulated regressor coefficient values in the order they appear in the
#'     latent state transition equation from d=1,...,DD
#'   }
#' @param par_inits list of initialization values for parameter s
#' @param par_names list of 3: first component is the name for sigma_squared and
#'   phi, the second component is a list of regressor names for Z-type
#'   regressors, and the last is for U-type regressors
#' @param par_trues if a simulation study, then one can provide true parameter
#'   values; the input argument should have the same structure as the argument
#'   \code{par_inits}
#' @param state_names character vector of dimension NN giving the names of
#'   the US states
#' @param energy_type_names character vector of dimension DD giving the names of
#'   the energy types
#' @param NN number of cross sectional units
#' @param TT number of time periods
#' @param DD dimension of the measurements (latent states)
#'
#' @return a list of 6 elements:
#' \itemize{
#'     \item{\code{$mcmc_sims:}}{matrix of dimension DD x num_mncmc_draws} of
#'     simulated standard deviation values
#'     \item{\code{$states:}}{matrix of dimension DD x num_mncmc_draws} of
#'     simulated autoregressive parameter values
#'     \item{\code{$par_names}}{matrix of dimension \code{num_par_beta x}
#'     \code{ num_mncmc_draws} of simulated regressor coefficient values in the
#'     order they appear in the latent state transition equation from
#'     d=1,...,DD}
#'   }
#'
#' @export
pgas_out_2_diagnostics <- function(pgas_out, par_inits,
                                   par_names = NULL, par_trues = NULL,
                                   state_names = NULL, energy_type_names = NULL,
                                   NN, TT, DD) {
  if (is.null(state_names)) state_names <- as.character(1:NN)
  if (is.null(energy_type_names)) energy_type_names <- as.character(1:DD)

  # browser()
  bet_z_dim <- sapply(par_inits[[3]][1:DD], length)
  if (!is.null(par_inits[4][[1]])) {
    bet_u_dim <- sapply(par_inits[[4]][1:DD], nrow)
    Unull <- FALSE
  } else {
    bet_u_dim <- 0
    Unull <- TRUE
  }
  if (is.null(par_names)) {
    par_names_sig_phi <- vector("list", 2)
    par_names_sig_phi[[1]] <- paste0(names(pgas_out)[1], "_", energy_type_names)
    par_names_sig_phi[[2]] <- paste0(names(pgas_out)[2], "_", energy_type_names)
    par_names_z <- list(DD)
    if (!Unull) par_names_u <- rep(list(vector("list", DD)), times = NN)
    for (d in 1:DD) {
      par_names_z[[d]] <- paste0("bet_z", 1:bet_z_dim[d], "_d=", d)
    }
    if (!Unull) {
      for (n in 1:NN) {
        for (d in 1:DD) {
          par_names_u[[n]][[d]] <- paste0("bet_u", 1:bet_u_dim[d],
                                          "_d=",  rep(d, times = bet_u_dim[d]),
                                          "_n=", n)
        }
      }
    }
  } else if (!is.null(par_names) && length(par_names) == 3) {
    par_names_sig_phi <- list(paste0(par_names[["reg_names_sig_phi"]][1], "_",
                                     energy_type_names),
                              paste0(par_names[["reg_names_sig_phi"]][2], "_",
                                     energy_type_names))
    par_names_z <- list(DD)
    if(!Unull) par_names_u <- rep(list(vector("list", DD)), times = NN)
    for (d in 1:DD) {
      par_names_z[[d]] <- paste0(energy_type_names[d], "_",
                                 par_names[["reg_names_Z"]][[d]])
    }
    if (!Unull) {
      for (n in 1:NN) {
        for (d in 1:DD) {
          for (i in 1:bet_u_dim[d]) {
            par_names_u[[n]][[d]] <- paste0(energy_type_names[d],
                                            "_", state_names[n],
                                            "_",
                                            par_names[["reg_names_U"]][[d]][[i]])
          }
        }
      }
    }
  } else {
    stop(paste0("Assignment of parameter names, ",
                "labels etc is not working; check internal code."))
  }
  # browser()
  num_pars  <- 2*DD + sum(bet_z_dim) + sum(bet_u_dim) * NN
  out        <- vector("list", 8)
  names(out) <- c("mcmc_sims",
                  "states",
                  "par_names",
                  "par_names_plots",
                  "lab_names",
                  "start_vals",
                  "true_vals")

  par_names_all       <- character(num_pars)
  lab_names_all       <- character(num_pars)
  par_names_all_plots <- character(num_pars)

  par_names_all <- c(par_names_sig_phi[[1]],  par_names_sig_phi[[2]],
                     unlist(par_names_z))
  if (!Unull) {
    par_names_all <- c(par_names_all, unlist(par_names_u))
  }
  lab_names_all <- par_names_all
  par_names_all_plots <- par_names_all
  # for (i in 1:2) {
  #   par_names_all[1:DD + DD*(i - 1)]       <- par_names_sig_phi[[i]]
  #   lab_names_all[1:DD + DD*(i - 1)]       <- par_names_sig_phi[[i]]
  #   par_names_all_plots[1:DD + DD*(i - 1)] <- par_names_sig_phi[[i]]
  # }
  # for (d in 1:DD) {
  #   id_start <- 2*DD + sum(c(1, bet_z_dim)[1:d])
  #   id_end <- 2*DD + sum(bet_z_dim[1:d])
  #   par_names_all[id_start:id_end]       <- par_names_z[[d]]
  #   lab_names_all[id_start:id_end]       <- par_names_z[[d]]
  #   par_names_all_plots[id_start:id_end] <- par_names_z[[d]]
  # }
  # if (!Unull) {
  #   browser()
  #   for (n in 1:NN) {
  #     for (d in 1:DD) {
  #       id_start <- 2*DD + sum(bet_z_dim) + sum(c(1, bet_u_dim)[1:d]) + sum(bet_u_dim[1:DD])*(n - 1)
  #       id_end   <- 2*DD + sum(bet_z_dim) + sum(bet_u_dim[1:d]) + sum(bet_u_dim[1:DD])*(n - 1)
  #       par_names_all[id_start:id_end]       <- par_names_u[[d]][[n]][d]
  #       lab_names_all[id_start:id_end]       <- par_names_u[[d]][[n]][d]
  #       par_names_all_plots[id_start:id_end] <- par_names_u[[d]][[n]][d]
  #     }
  #   }
  # }
  out$num_pars  <- num_pars
  mcmc_sims_tmp <- cbind(t(pgas_out[[1]]), t(pgas_out[[2]]), t(pgas_out[[3]]))
  if (!Unull) {
    for (n in 1:NN) {
      # browser()
      mcmc_sims_tmp <- cbind(mcmc_sims_tmp, t(pgas_out[[4]][, , n]))
    }
  }
  out$mcmc_sims       <- mcmc_sims_tmp
  out$states          <- pgas_out[["x"]]
  out$par_names       <- par_names_all
  out$par_names_plots <- par_names_all_plots
  out$lab_names       <- lab_names_all
  if (!Unull) {
    out$start_vals      <- unname(c(unlist(par_inits[1:3]),
                                    as.vector(Reduce(function(x, y) {
                                      return(rbind(x,y))},
                                      par_inits[[4]]))))
  } else {
    out$start_vals      <- unname(c(unlist(par_inits[1:3])))
  }
  if (!is.null(par_trues)) {
    out$true_vals  <- c(par_trues$sig_sq[, 1, drop = TRUE],
                        par_trues$phi[, 1, drop = TRUE],
                        unlist(par_trues$bet_z))
    if (!Unull) {
      out$true_vals <- c(out$true_vals, as.vector(
        Reduce(function(x, y) {return(rbind(x,y))},
               par_trues[[4]])))
    }
  }
  return(out)
}
#' Transforms PGAS output to list format for comparison/testing reasons
#'
#' @param pgas_out output of \code{pgas_cpp()} or \code{pgas_R()}
#' @param NN cross sectional dimension
#' @param TT number of time periods
#' @param DD multivariate response/measurement dimension: e.g. number of
#'   shares/fractions if the measurements are from a Dirichlet
#' @param MM number of overall (particle) MCMC iterations
#' @param dim_bet_z number of regressors coefficients of z-type regressors for
#'   each \code{d,...,DD}
#' @param cpp logical; if \code{TRUE}, the particle gibbs output comes from cpp
#'   code, if \code{FALSE} from R code
#' @return a list of 19 elements corresponding to former output format
#' @export
pgas_out_2_list <- function(pgas_out, DD, NN, TT, MM, dim_bet_z, cpp = NULL) {
  # browser()
  if (is.null(cpp)) {
    stop("Specify if pgas output comes from a pgas_cpp() or from pgas_R()!")
  }
  id_bet <- c(1, cumsum(dim_bet_z) + 1)
  out   <- list()
  xtraj <- list()
  for (n in 1:NN) {
    for (d in 1:DD) {
      out[[1 + (d - 1)*3]] <- pgas_out[[1]][d, ]
      out[[2 + (d - 1)*3]] <- pgas_out[[2]][d, ]
      out[[3 + (d - 1)*3]] <- pgas_out[[3]][id_bet[d]:(id_bet[d + 1] - 1), ]
      if (cpp) {
        xtraj[[d]] <- pgas_out[[4]][(1 + TT*(d - 1)):(TT*d), , n]
      } else {
        xtraj[[d]] <- t(pgas_out[[4]][ , d, , n])
      }
    }
  }
  out[[19]] <- xtraj

  names(out) <-  paste0(c("sigma_sq_xa",
                          "phi_xa",
                          "bet_xa"), 1:DD)
  return(out)
}
