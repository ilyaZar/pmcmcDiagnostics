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
#' @param mcmc_sims simulated draws/ mcmc output: matrix of dimension "simulated
#'   MCMC draws" (rows) x "parameters" (cols)
#' @param states state trajectories i.e. particle filter (SMC) output
#' @param par_names parameter names as used in the pgas output (col names for
#'   \code{mcmcm_sims} e.g.)
#' @param par_names_plots parameter names for plot labelling
#' @param lab_names names for labels
#' @param start_vals start values of the parameters
#' @param true_vals true values if a simulation is run; default is NULL and true
#'   values will not be added for the histogram and trace plots
#' @param burn burnin period
#' @param thin thinning; \code{thin = 10} means every 10th draw is excluded
#' @param ki_prob probablitlity with which confidence intervals and highest
#'   posterior density regions are computed
#' @param plot_view logical; if \code{TRUE}, plots are returned for viewing
#'   (RStudio pane)
#' @param plot_ggp2 logical; if \code{TRUE},ggplots are generated (if
#'   \code{FALSE}, then base plots are used)
#' @param plot_save logical; if \code{TRUE}, plots are saved with names
#'   specified under \code{plot_name} and path given by the argument
#'   \code{plot_path}
#' @param plot_name if \code{plot_save = TRUE}, defines the names of the plots
#'   as a character vector
#' @param plot_path if \code{plot_save = TRUE}, defines the path where plots are
#'   stored as a character
#' @param table_view logical; if \code{TRUE}, output table can be viewd (RStudio
#'   pane)
#' @param table_save logical; if \code{TRUE}, output tables are saved
#' @param table_name if \code{table_view = TRUE}, defines the name of the output
#'   table as a character
#' @param table_path if \code{table_view = TRUE}, defines the path where output
#'   table is stored as character vector
#' @param table_prec if \code{table_view = TRUE}, rounding digits for values of
#'   the output table
#' @param compute_ess logical; if \code{TRUE}, computes the effective sample
#'   size
#' @param compute_ess_stan logical; if \code{TRUE}, computes the effective
#'   sample sizes in terms of ESS bulk and ESS tail (see the function
#'   documentation of \code{diagnostics_table()} for details)
#' @param ur_view logical; if \code{TRUE} are returned for viewing (RStuido
#'   pane)
#' @param ur_save logical; if \code{TRUE} update rates are saved
#' @param ur_name if \code{ur_save == TRUE}, specifies name of the update rate
#'   plot
#' @param ur_path if \code{ur_save == TRUE}, specifies path of the update rate
#'   plot
#'
#' @return nothing but plots and table are saved and displayed upon request
#' @export
analyse_mcmc_convergence2 <- function(mcmc_sims, states,
                                      par_names, par_names_plots,
                                      lab_names = NULL,
                                      start_vals, true_vals = NULL,
                                      burn, thin = NULL,
                                      ki_prob   = 0.9,
                                      plot_view = FALSE,
                                      plot_ggp2 = FALSE,
                                      plot_save = FALSE,
                                      plot_name = "",
                                      plot_path = NULL,
                                      table_view = FALSE,
                                      table_save = FALSE,
                                      table_name = "",
                                      table_path = NULL,
                                      table_prec = 4,
                                      compute_ess = TRUE,
                                      compute_ess_stan = FALSE,
                                      ur_view = FALSE,
                                      ur_save = FALSE,
                                      ur_name = "",
                                      ur_path = NULL) {
  num_mcmc        <- dim(mcmc_sims)[1]
  num_par         <- dim(mcmc_sims)[2]

  if (num_mcmc - burn < 1) {
    stop("Burn-in period too large: the number of (P)MCMC samples is: ",
         num_mcmc, " while burn-in is: ", burn, "!")
  }

  mcmc_sims_after <- burn_and_thin(mcmc_sims, burnin = burn, thin)

  mcmc_sims_df        <- data.frame(cbind(1:num_mcmc, mcmc_sims))
  names(mcmc_sims_df) <- c("num_mcmc", par_names)
  mcmc_sims_df_after  <- subset(mcmc_sims_df, num_mcmc >= burn)

  posterior_means <- colMeans(mcmc_sims_after)
  #
  #
  #
  #
  #
  if (plot_view || plot_save) {
    browser()
    for (i in 1:num_par) {
      if (plot_ggp2) {
        plot_returned <- generate_ggplot2(mcmc_sims_df = mcmc_sims_df,
                                          mcmc_sims_df_after=mcmc_sims_df_after,
                                          burn = burn,
                                          thin = thin,
                                          num_mcmc = num_mcmc,
                                          par_names = par_names,
                                          true_vals = true_vals,
                                          posterior_means = posterior_means,
                                          plot_num = i)
        if (plot_view) {
          plot_current <- gridExtra::grid.arrange(plot_returned[[1]],
                                                  plot_returned[[2]],
                                                  plot_returned[[3]],
                                                  plot_returned[[4]],
                                                  nrow = 2,
                                                  top = grid::textGrob(par_names_plots[i],
                                                           gp = grid::gpar(fontsize = 20, font = 3)))
          print(plot_current)
        }
        if (plot_save) {
          plot_current <- gridExtra::arrangeGrob(plot_returned[[1]],
                                                 plot_returned[[2]],
                                                 plot_returned[[3]],
                                                 plot_returned[[4]],
                                                 nrow = 2,
                                                 top = grid::textGrob(par_names_plots[i],
                                                                      gp = grid::gpar(fontsize = 20, font = 3)))
          ggplot2::ggsave(filename = paste(plot_name, "_", par_names[i], ".pdf", sep = ""),
                          plot = plot_current,
                          path = plot_path,
                          device = "eps",
                          width = 18,
                          height = 10.5,
                          units = "cm")
          current_plot_name <- file.path(plot_path,
                                         paste(plot_name,"_",
                                               par_names[i],
                                               ".pdf",
                                               sep = ""))
          print(paste("Saved plots in: ", current_plot_name))
        }
      } else {
        if (plot_view) {
          plot_returned <- generate_plot2(mcmc_sims = mcmc_sims,
                                          mcmc_sims_after = mcmc_sims_after,
                                          burn = burn,
                                          thin = thin,
                                          num_mcmc = num_mcmc,
                                          par_names = par_names,
                                          par_names_plots = par_names_plots,
                                          true_vals = true_vals,
                                          posterior_means = posterior_means,
                                          plot_num = i)
        }
        if (plot_save) {
          current_plot_name <- file.path(plot_path,
                                         paste(plot_name, "_",
                                               par_names[i],
                                               ".pdf",
                                               sep = ""))
          grDevices::setEPS()
          grDevices::postscript(current_plot_name, width = 18, height = 10.5)
          generate_plot2(mcmc_sims = mcmc_sims,
                         mcmc_sims_after = mcmc_sims_after,
                         burn = burn,
                         num_mcmc = num_mcmc,
                         par_names = par_names,
                         par_names_plots = par_names_plots,
                         true_vals = true_vals,
                         posterior_means = posterior_means,
                         plot_num = i)
          grDevices::dev.off()
          print(paste("Saved plots in: ", current_plot_name))
        }
      }
    }
  }
  #
  #
  #
  #
  #
  browser()
  summary_results <- diagnostics_table(num_par = num_par,
                                       mcmc_sims = mcmc_sims,
                                       mcmc_sims_after = mcmc_sims_after,
                                       burn = burn,
                                       num_mcmc = num_mcmc,
                                       posterior_means = posterior_means,
                                       start_vals = start_vals,
                                       true_vals = true_vals,
                                       ki_prob = ki_prob,
                                       compute_ess = compute_ess,
                                       compute_ess_stan = compute_ess_stan)
  if (table_view) {
    summary_results_view <- summary_results
    ID_round <- which(!sapply(summary_results_view, is.logical))
    summary_results_view[, ID_round] <- round(summary_results_view[, ID_round],
                                              digits = table_prec)
    row.names(summary_results_view) <- par_names
    utils::View(summary_results_view, title = paste(table_name,
                                             "_summary_results",
                                             sep = ""))
  }
  if (table_save) {
    if (is.null(par_names_plots) || is.null(lab_names)) {
      stop("Can't save results in table form: label names required!")
    }
    summary_results_save <- summary_results
    summary_results_save <- cbind(lab_names, summary_results_save)
    # row.names(summary_results_save) <- par_names
    # summary_results_save <- cbind(par_name = par_names, summary_results_save)
    readr::write_csv(summary_results_save, path = file.path(table_path,
                                                            paste0(table_name,
                                                                   ".csv")))
  }
  #
  #
  #
  #
  #
  if (ur_view) {
    graphics::par(mfrow = c(1, 1))
    analyse_states_ur(trajectories = states)
  }
  if (ur_save) {
    current_plot_name <- file.path(plot_path, paste0("00_", ur_name, ".pdf"))
    grDevices::setEPS()
    grDevices::postscript(current_plot_name, width = 9, height = 5.25)
    analyse_states_ur(trajectories = states)
    grDevices::dev.off()
    print(paste("Saved update rate plot in: ", current_plot_name))
  }
  if (table_view) {
    return(summary_results_view)
  } else {
    return(invisible(summary_results_view))
  }
}
