#' Displays output diagnostics
#'
#' @description
#' Displays output diagnostic plots and tables of (particle) MCMC parameter
#' draws: histogram, trace and autocorrelations plots for all parameters as
#' (ggplot2 or base) plots, and an output table summary with posterior mean,
#' standard deviation of the parameter under the posterior distribution and
#' standard deviation of the posterior mean, confidence bands, HPDs, and
#' effective sample sizes (the latter also in the Stan variant i.e. ESS bulk and
#' tail).
#'
#' Moreover, update rates of simulated latent state processes (from e.g. a
#' bootstrap particle filter or general SMC output) are returned.
#'
#' @details
#' All diagnostics (plots/tables) can be saved when specifying a directory name
#' and setting the corresponding logical flag \code{XXX_save = TRUE} (where
#' \code{XXX} can be \code{plot}, \code{table} etc., indicating whether plots,
#' output diagnostics tables or update rates should be saved).
#'
#' @param model_output a model as produced from the output of [BNMPD::pgas_d]
#'   e.g.
#' @param mcmc_sims simulated draws/MCMC output: matrix of dimension "simulated
#'   MCMC draws" (rows) x "parameters" (cols)
#' @param states state trajectories i.e. particle filter (SMC) output
#' @param model_meta a list of two elements:
#'   \itemize{
#'     \item{\code{par_val_names = NULL: }}{parameter names as used in the PGAS
#'     output (col names for \code{mcmcm_sims} e.g.)}
#'     \item{\code{par_lab_names = NULL: }}{names for labels}
#'   }
#' @param settings_mcmc a list of four elements:
#'   \itemize{
#'     \item{\code{burn: }}{burnin period - defaults to half of the MCMC draws}
#'     \item{\code{thin: }}{thinning; \code{thin = 10} means every 10th draw is
#'     excluded}
#'     \item{\code{ki_prob: }}{probability with which confidence intervals and
#'     highest posterior density regions are computed}
#'     \item{\code{compute_ess: }}{logical; if \code{TRUE}, computes the effective sample
#'   size}
#'   \item{\code{compute_ess_stan: }}{logical; if \code{TRUE}, computes the
#'   effective sample sizes in terms of ESS bulk and ESS tail (see the function
#'   documentation of \code{diagnostics_table()} for details)}
#'   }
#' @param settings_plots a list of five:
#'   \itemize{
#'     \item{\code{plot_view: }}{logical; if \code{TRUE}, plots are returned for
#'     viewing (RStudio pane)}
#'     \item{\code{plot_ggp2: }}{logical; if \code{TRUE},ggplots are generated (if
#'     \code{FALSE}, then base plots are used)}
#'     \item{\code{plot_save: }}{logical; if \code{TRUE}, plots are saved with
#'     names equal to parameter name and path given via \code{plot_path}}
#'     \item{\code{plot_save_all: }}{logical; if \code{TRUE}, plots are saved with
#'     name specified under \code{plot_name} and path given given via
#'     \code{plot_path}}
#'     \item{\code{plot_name: }}{if \code{plot_save = TRUE}, defines the names of
#'     the plots as a character vector}
#'     \item{\code{plot_path: }}{if \code{plot_save = TRUE}, defines the path
#'     where plots are stored as a character}
#'   }
#' @param settings_table a list of five:
#'   \itemize{
#'     \item{\code{table_view: }}{logical; if \code{TRUE}, output table can be
#'     viewed (RStudio pane)}
#'     \item{\code{table_save: }}{logical; if \code{TRUE}, output tables are
#'     saved}
#'     \item{\code{table_name: }}{if \code{table_view = TRUE}, defines the name of
#'     the output table as a character}
#'     \item{\code{table_path: }}{ if \code{table_view = TRUE}, defines the path
#'     where output table is stored as character vector}
#'     \item{\code{table_prec: }}{if \code{table_view = TRUE}, rounding digits for
#'     values of the output table}
#'   }
#' @param settings_urs a list of four:
#'   \itemize{
#'     \item{\code{ur_view: }}{logical; if \code{TRUE} are returned for viewing
#'     (RStudio pane)}
#'     \item{\code{ur_save: }}{logical; if \code{TRUE} update rates are saved}
#'     \item{\code{ur_name: }}{if \code{ur_save == TRUE}, specifies name of the
#'     update rate plot}
#'     \item{\code{ur_path: }}{if \code{ur_save == TRUE}, specifies path of the
#'     update rate plot}
#'   }
#' @return nothing but plots and table are saved and displayed upon request
#' @export
analyse_mcmc_convergence2 <- function(model_output = NULL,
                                      mcmc_sims = NULL,
                                      states = NULL,
                                      model_meta = list(par_val_names = NULL,
                                                        par_lab_names = NULL),
                                      settings_mcmc = list(burn = round(nrow(mcmc_sims)/2,
                                                                        digits = 0),
                                                           thin = NULL,
                                                           ki_prob = 0.9,
                                                           q_probs = c(0.025,
                                                                       0.25,
                                                                       0.5,
                                                                       0.75,
                                                                       0.975),
                                                           compute_ess = TRUE,
                                                           compute_ess_stan = TRUE),
                                      settings_plots = list(plot_view = FALSE,
                                                            plot_ggp2 = FALSE,
                                                            plot_save = FALSE,
                                                            plot_save_all = FALSE,
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
                  (!is.null(mcmc_sims) &&
                     !is.null(states))))

  if (inherits(x = model_output, what = "pmcmc")) states <- model_output$x

  par_names <- unname(unlist(model_meta$par_val_names))
  lab_names <- unname(unlist(model_meta$par_lab_names))
  true_vals <- get_true_vals(list_true_vals = model_output$true_vals)

  mcmc_sims <- model_out2sims(model_output, lab_names)
  num_mcmc <- dim(mcmc_sims)[1]
  num_par  <- dim(mcmc_sims)[2]

  if (num_mcmc - settings_mcmc$burn < 1) {
    stop("Burn-in period too large: the number of (P)MCMC samples is: ",
         num_mcmc, " while burn-in is: ",
         settings_mcmc$burn, "!")
  }

  start_vals         <- mcmc_sims[1, ]
  mcmc_sims_after    <- burn_and_thin(mcmc_sims,
                                      burnin = settings_mcmc$burn,
                                      settings_mcmc$thin)
  mcmc_sims_df       <- data.frame(cbind(num_mcmc = 1:num_mcmc, mcmc_sims))
  mcmc_sims_df_after <- subset(mcmc_sims_df,
                               num_mcmc >= settings_mcmc$burn)
  posterior_means    <- colMeans(mcmc_sims_after)

  if (settings_plots$plot_view ||
      settings_plots$plot_save ||
      settings_plots$plot_save_all) {
    if (settings_plots$plot_ggp2) {
      generate_ggplot2_all(mcmc_sims_df, mcmc_sims_df_after,
                           settings_mcmc$burn,
                           settings_mcmc$thin,
                           num_mcmc,
                           lab_names, true_vals,
                           posterior_means,
                           settings_plots,
                           lab_names)
    } else {
      generate_plot_all(mcmc_sims, mcmc_sims_after,
                        settings_mcmc$burn, settings_mcmc$thin,
                        num_mcmc,
                        lab_names, true_vals,
                        posterior_means,
                        settings_plots,
                        lab_names)
    }
  }
  out_mcmc <- NULL
  if (settings_table$table_view || settings_table$table_save) {
    out_mcmc <- diagnostics_table(num_par = num_par,
                                  par_names = lab_names,
                                  mcmc_sims = mcmc_sims,
                                  mcmc_sims_after = mcmc_sims_after,
                                  burn = settings_mcmc$burn,
                                  num_mcmc = num_mcmc,
                                  posterior_means = posterior_means,
                                  start_vals  = start_vals ,
                                  true_vals = true_vals,
                                  ki_prob = settings_mcmc$ki_prob,
                                  q_probs = settings_mcmc$q_probs,
                                  ESS_STANDARD = settings_mcmc$compute_ess,
                                  ESS_STAN = settings_mcmc$compute_ess_stan,
                                  settings_table)
  }

  out_urs <- NULL
  if (settings_urs$ur_view || settings_urs$ur_save) {
    out_urs <- get_update_rates(states, settings_urs, settings_plots)
  }

  if (is.null(out_mcmc)) warning(paste0("No MCMC summary results computed, ",
                                        "nothing to return ..."))
  if (is.null(out_urs)) warning(paste0("No update rates computed, ",
                                       "nothing to return ..."))
  return(list(out_mcmc = out_mcmc,
              out_urs = out_urs))
}
