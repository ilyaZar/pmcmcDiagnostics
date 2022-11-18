#' Dispays output diagnostics
#'
#' Displays output diagnostic plots and tables of (particle) MCMC parameter
#' draws: histogram, trace and autocorrelations plots for all parameters as
#' (ggplot2 or base) plots, and an output table summary with posterior mean,
#' standard deviation of the parameter under the posterior distribution and
#' standard deviation of the posterior mean, confidence bands, HPDs, and
#' effective sample sizes (the latter also in the stan variant i.e. ESS bulk and
#' tail). Also, update rates of simulated latent state parameters (bootsrap
#' particle filter/SMC output) are returned. All diagnostics can be saved
#' specifiying a directory, name and if a logical value, \code{..._save},
#' idnicating whether plots, output diagnostics table or update rates should be
#' save.
#'
#' @param model_output a model as produced from the output of [BNMPD::pgas_d]
#'   e.g.
#' @param mcmc_sims simulated draws/ mcmc output: matrix of dimension "simulated
#'   MCMC draws" (rows) x "parameters" (cols)
#' @param states state trajectories i.e. particle filter (SMC) output
#' @param model_meta a list of two elements:
#'   \itemize{
#'     \item{\code{par_val_names = NULL}}{parameter names as used in the pgas
#'     output (col names for \code{mcmcm_sims} e.g.)}
#'     \item{\code{par_lab_names = NULL}}{names for labels}
#'   }
#' @param settings_mcmc a list of four elements:
#'   \itemize{
#'     \item{\code{burn}}{burnin period}
#'     \item{\code{thin}}{thinning; \code{thin = 10} means every 10th draw is
#'     excluded}
#'     \item{\code{ki_prob}}{ probablitlity with which confidence intervals and
#'     highest posterior density regions are computed}
#'     \item{\code{compute_ess}}{logical; if \code{TRUE}, computes the effective sample
#'   size}
#'     \item{\code{compute_ess}}{logical; if \code{TRUE}, computes the effective
#'   sample sizes in terms of ESS bulk and ESS tail (see the function
#'   documentation of \code{diagnostics_table()} for details)}
#'   }
#' @param settings_plots a list of five:
#'   \itemize{
#'     \item{\code{plot_view}}{logical; if \code{TRUE}, plots are returned for
#'     viewing (RStudio pane)}
#'     \item{\code{plot_ggp2}}{logical; if \code{TRUE},ggplots are generated (if
#'   \code{FALSE}, then base plots are used)}
#'     \item{\code{pplot_save}}{logical; if \code{TRUE}, plots are saved with
#'     names specified under \code{plot_name} and path given by the argument
#'   \code{plot_path}}
#'     \item{\code{plot_name}}{if \code{plot_save = TRUE}, defines the names of
#'     the plots as a character vector}
#'     \item{\code{plot_path}}{if \code{plot_save = TRUE}, defines the path
#'     where plots are stored as a character}
#'   }
#' @param settings_table a list of five:
#'   \itemize{
#'     \item{\code{table_view}}{ logical; if \code{TRUE}, output table can be
#'     viewed (RStudio pane)}
#'     \item{\code{table_save}}{logical; if \code{TRUE}, output tables are
#'     saved}
#'     \item{\code{table_name}}{if \code{table_view = TRUE}, defines the name of
#'     the output table as a character}
#'     \item{\code{table_path}}{ if \code{table_view = TRUE}, defines the path
#'     where output table is stored as character vector}
#'     \item{\code{table_prec}}{if \code{table_view = TRUE}, rounding digits for
#'     values of the output table}
#'   }
#' @param settings_urs a list of four:
#'   \itemize{
#'     \item{\code{ur_view}}{logical; if \code{TRUE} are returned for viewing
#'     (RStudio pane)}
#'     \item{\code{ur_save}}{logical; if \code{TRUE} update rates are saved}
#'     \item{\code{ur_name}}{if \code{ur_save == TRUE}, specifies name of the
#'     update rate plot}
#'     \item{\code{ur_path}}{if \code{ur_save == TRUE}, specifies path of the
#'     update rate plot}
#'   }
#' @return nothing but plots and table are saved and displayed upon request
#' @export
analyse_mcmc_convergence2 <- function(model_output = NULL,
                                      mcmc_sims = NULL,
                                      states = NULL,
                                      model_meta = list(par_val_names = NULL,
                                                        par_lab_names = NULL),
                                      settings_mcmc = list(burn = 1,
                                                           thin = NULL,
                                                           ki_prob = 0.9,
                                                           compute_ess = TRUE,
                                                           compute_ess_stan = FALSE),
                                      settings_plots = list(plot_view = FALSE,
                                                            plot_ggp2 = FALSE,
                                                            plot_save = FALSE,
                                                            plot_name = "",
                                                            plot_path = NULL),
                                      settings_table = list(table_view = FALSE,
                                                            table_save = FALSE,
                                                            table_name = "",
                                                            table_path = NULL,
                                                            table_prec = 4),
                                      settings_urs = list(ur_view = FALSE,
                                                          ur_save = FALSE,
                                                          ur_name = "",
                                                          ur_path = NULL)) {
  stopifnot(any(!is.null(model_output) ||
                  !is.null(mcmc_sims) ||
                  !is.null(states)))
  summary_results_view <- NULL
  if (class(model_output) == "pmcmc") {
    states <- model_output$x
  }
  par_names <- unlist(model_meta$par_val_names)
  lab_names <- unlist(model_meta$par_lab_names)
  mcmc_sims <- model_out2sims(model_output, par_names)
  true_vals <- get_true_vals(list_true_vals = model_output$true_vals)
  # par_names,
  # lab_names = NULL,
  num_mcmc <- dim(mcmc_sims)[1]
  num_par  <- dim(mcmc_sims)[2] - 1
  burn <- settings_mcmc$burn
  thin <- settings_mcmc$thin

  if (num_mcmc - burn < 1) {
    stop("Burn-in period too large: the number of (P)MCMC samples is: ",
         num_mcmc, " while burn-in is: ", burn, "!")
  }
  start_vals         <- mcmc_sims[1, ]
  mcmc_sims_after    <- burn_and_thin(mcmc_sims, burnin = burn, thin)
  mcmc_sims_df       <- data.frame(cbind(num_mcmc = 1:num_mcmc, mcmc_sims))
  mcmc_sims_df_after <- subset(mcmc_sims_df, num_mcmc >= burn)
  posterior_means    <- colMeans(mcmc_sims_after[,-1])
  #
  #
  #
  #
  #
  if (settings_plots$plot_view || settings_plots$plot_save) {
    for (i in 1:num_par) {
      if (settings_plots$plot_ggp2) {
        plot_returned <- generate_ggplot2(mcmc_sims_df = mcmc_sims_df[, -1],
                                          mcmc_sims_df_after=mcmc_sims_df_after[, -1],
                                          burn = burn,
                                          thin = thin,
                                          num_mcmc = num_mcmc,
                                          par_names = par_names,
                                          true_vals = true_vals,
                                          posterior_means = posterior_means,
                                          plot_num = i)
        if (settings_plots$plot_view) {
          plot_current <- gridExtra::grid.arrange(plot_returned[[1]],
                                                  plot_returned[[2]],
                                                  plot_returned[[3]],
                                                  plot_returned[[4]],
                                                  nrow = 2,
                                                  top = grid::textGrob(lab_names[i],
                                                           gp = grid::gpar(fontsize = 20, font = 3)))
          print(plot_current)
        }
        if (settings_plots$plot_save) {
          plot_current <- gridExtra::arrangeGrob(plot_returned[[1]],
                                                 plot_returned[[2]],
                                                 plot_returned[[3]],
                                                 plot_returned[[4]],
                                                 nrow = 2,
                                                 top = grid::textGrob(lab_names[i],
                                                                      gp = grid::gpar(fontsize = 20, font = 3)))
          ggplot2::ggsave(filename = paste(settings_plots$plot_name, "_", par_names[i], ".pdf", sep = ""),
                          plot = plot_current,
                          path = settings_plots$plot_path,
                          device = "eps",
                          width = 18,
                          height = 10.5,
                          units = "cm")
          current_plot_name <- file.path(settings_plots$plot_path,
                                         paste(settings_plots$plot_name,"_",
                                               par_names[i],
                                               ".pdf",
                                               sep = ""))
          print(paste("Saved plots in: ", current_plot_name))
        }
      } else {
        if (settings_plots$plot_view) {
          plot_returned <- generate_plot2(mcmc_sims = mcmc_sims[, -1],
                                          mcmc_sims_after = mcmc_sims_after[, -1],
                                          burn = burn,
                                          thin = thin,
                                          num_mcmc = num_mcmc,
                                          par_names = par_names,
                                          par_names_plot = lab_names,
                                          true_vals = true_vals,
                                          posterior_means = posterior_means,
                                          plot_num = i)
        }
        if (settings_plots$plot_save) {
          current_plot_name <- file.path(settings_plots$plot_path,
                                         paste(settings_plots$plot_name, "_",
                                               par_names[i],
                                               ".pdf",
                                               sep = ""))
          grDevices::setEPS()
          grDevices::postscript(current_plot_name, width = 18, height = 10.5)
          generate_plot2(mcmc_sims = mcmc_sims[, -1],
                         mcmc_sims_after = mcmc_sims_after[, -1],
                         burn = burn,
                         num_mcmc = num_mcmc,
                         par_names = par_names,
                         par_names_plot = lab_names,
                         true_vals = true_vals,
                         posterior_means = posterior_means,
                         plot_num = i)
          grDevices::dev.off()
          print(paste("Saved plots in: ", current_plot_name))
        }
      }
    }
  }
  browser()
  summary_results <- diagnostics_table(num_par = num_par,
                                       mcmc_sims = mcmc_sims[, -1],
                                       mcmc_sims_after = mcmc_sims_after[, -1],
                                       burn = burn,
                                       num_mcmc = num_mcmc,
                                       posterior_means = posterior_means,
                                       start_vals  = start_vals ,
                                       true_vals = true_vals,
                                       ki_prob = settings_mcmc$ki_prob,
                                       compute_ess = settings_mcmc$compute_ess,
                                       compute_ess_stan = settings_mcmc$compute_ess_stan)
  if (settings_table$table_view) {
    summary_results_view <- summary_results
    ID_round <- which(!sapply(summary_results_view, is.logical))
    summary_results_view[, ID_round] <- round(summary_results_view[, ID_round],
                                              digits = settings_table$table_prec)
    row.names(summary_results_view) <- par_names
    utils::View(summary_results_view, title = paste(settings_table$table_name,
                                             "_summary_results",
                                             sep = ""))
  }
  if (settings_table$table_save) {
    if (is.null(lab_names) || is.null(lab_names)) {
      stop("Can't save results in table form: label names required!")
    }
    summary_results_save <- summary_results
    summary_results_save <- cbind(lab_names, summary_results_save)
    # row.names(summary_results_save) <- par_names
    # summary_results_save <- cbind(par_name = par_names, summary_results_save)
    readr::write_csv(summary_results_save,
                     file = file.path(settings_table$table_path,
                                      paste0(settings_table$table_name,
                                             ".csv")))
  }
  #
  #
  #
  #
  #
  if (settings_urs$ur_view) {
    graphics::par(mfrow = c(1, 1))
    analyse_states_ur(trajectories = states)
  }
  if (settings_urs$ur_save) {
    current_plot_name <- file.path(settings_plots$plot_path,
                                   paste0("00_", settings_urs$ur_name, ".pdf"))
    grDevices::setEPS()
    grDevices::postscript(current_plot_name, width = 9, height = 5.25)
    analyse_states_ur(trajectories = states)
    grDevices::dev.off()
    print(paste("Saved update rate plot in: ", current_plot_name))
  }
  if (!is.null(summary_results_view)) {
    return(summary_results_view)
  } else {
    warning("No summary results computed, nothing to return ...")
    return(invisible(summary_results_view))
  }
}
model_out2sims <- function(mod_out, par_names) {
  sig_sq_x  <- mod_out$sig_sq_x
  phi_x     <- mod_out$phi_x
  bet_z     <- mod_out$bet_z
  bet_u     <- mod_out$bet_u
  vcm_bet_u <- mod_out$vcm_bet_u


  if (!is.null(sig_sq_x)) {
    sig_sq_x <- matrix(sig_sq_x, ncol = 1)
  }
  if (!is.null(phi_x)) {
    phi_x <- matrix(phi_x, ncol = 1)
  }
  if (!is.null(bet_z)) {
    if (!is.na(dim(bet_z)[3])) {
      bet_z <- matrix(apply(bet_z, 3, t), nrow = nrow(bet_z[,,1]))
    } else {
      bet_z <- t(bet_z)
    }
  }
  if (!is.null(bet_u)) {
    if (!is.na(dim(bet_u)[3])) {
      NN     <- dim(bet_u)[3]
      MM     <- dim(bet_u)[2]
      num_re <- dim(bet_u)[1]

      out <- matrix(0, ncol = num_re*NN, nrow = MM)
      iter <- 1
      for (r in 1:num_re) {
        for (n in 1:NN) {
          out[, iter] <- as.vector(bet_u[r, , n])
          iter <- iter + 1
        }
      }

     bet_u <- out
    } else {
      stop("not yet implemented")
    }
  }
  for (d in 1:length(mod_out$vcm_bet_u)) {
    dim_vcm_bet_u <- nrow(mod_out$vcm_bet_u[[d]])
    vcm_bet_u_mat <- matrix(0, nrow = dim(mod_out$vcm_bet_u[[1]])[3],
                            ncol = dim_vcm_bet_u^2)
    k <- 1
    for (i in 1:dim_vcm_bet_u) {
      for (j in 1:dim_vcm_bet_u) {
        vcm_bet_u_mat[, k] <- mod_out$vcm_bet_u[[1]][i, j, ]
        k <- k + 1
      }
    }
  }
  par_list <- list(sig_sq_x, phi_x, bet_z, bet_u, vcm_bet_u_mat)
  id_bind <- sapply(par_list, is.null)
  out <- Reduce(cbind,
                par_list[!id_bind],
                init =  1:nrow(par_list[!id_bind][[1]]))
  colnames(out) <- c("sim_mcmc", par_names)
  return(out)
}
get_true_vals <- function(list_true_vals) {
  names_pars <- names(list_true_vals)

  DD <- nrow(list_true_vals[["sig_sq"]])
  out <- numeric(0)
  for(j in names_pars) {
    for (d in 1:DD) {
      if (is.null(list_true_vals[[j]])) next;
      if (j %in% c("sig_sq", "phi")) {
        out <- c(out, unique(list_true_vals[[j]][d, ]))
      } else if (j == "beta_z_lin") {
        out <- c(out, unlist(list_true_vals[[j]][[d]]))
      } else if (j == "beta_u_lin") {
        out <- c(out, unlist(t(list_true_vals[[j]][[d]])))
      } else if (j == "vcm_u_lin") {
        out <- c(out, unlist(list_true_vals[[j]][[d]]))
      }
    }
  }
  return(out)
}
