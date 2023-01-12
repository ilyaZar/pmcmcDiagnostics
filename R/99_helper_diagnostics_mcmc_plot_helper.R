#' Generates base graphic analysis plot of MCMC output.
#'
#' Generates base graphic analysis plot of MCMC output for a model parameter
#' including: 1. histogram of mcmc draws (after burn-in and thinning) 2. trace
#' plot of mcmc draws (after burn-in and thinning) 3. autocorrelation plot of
#' mcmc draws (after burn-in and thinning) 4. trace plot of mcmc draws (BEFORE
#' burn-in and thinning)
#'
#' @param mcmc_sims mcmc draws as a matrix
#' @param mcmc_sims_after mcmc draws as a data.frame after burn-in and
#'   thinning
#' @param burn burn-in period
#' @param thin thinning; \code{thin = 10} means every 10th draw is excluded
#' @param num_mcmc total number of MCMC draws
#' @param par_names parameter names as used in the pgas output (col names for
#'   \code{mcmcm_sims} e.g.)
#' @param par_names_plots parameter names for plot labeling
#' @param true_vals true values if a simulation is run; default is NULL and true
#'   values will not be added for the histogram and trace plots
#' @param posterior_means posterior means (have to be pre-computed and passed
#'   directly here)
#' @param plot_num number of parameter or plot
#'
#' @return base plot in 2x2 form displaying the mcmc-diagnostics
generate_plot2 <- function(mcmc_sims,
                           mcmc_sims_after,
                           burn = 0,
                           thin = 0,
                           num_mcmc,
                           par_names,
                           par_names_plots,
                           true_vals = NULL,
                           posterior_means,
                           plot_num) {
  TRUE_VAL_AVAIL <- !is.null(true_vals) && !is.null(true_vals)

  graphics::par(mfrow = c(2, 2),
                oma = c(0,0,5,0))

  graphics::hist(mcmc_sims_after[, plot_num],
                 xlab = "values (after burn-in)",
                 main = "posterior density histogram")
  graphics::abline(v = posterior_means[plot_num], col = "red")
  if (TRUE_VAL_AVAIL) {
    graphics::abline(v = true_vals[plot_num], col = "green")
  }

  graphics::plot(mcmc_sims_after[, plot_num], type = "l",
                 xlab = "mcmc iterations (after burn-in)",
                 ylab = "sampled values",
                 main = paste("trace after burn-in", burn, sep = ": "))
  graphics::abline(h = posterior_means[plot_num], col = "red")
  if (TRUE_VAL_AVAIL) {
    graphics::abline(h = true_vals[plot_num], col = "green")
  }

  coda::autocorr.plot(mcmc_sims_after[, plot_num],
                      auto.layout = FALSE,
                      main = "autocorrelation")

  graphics::plot(mcmc_sims[, plot_num], type = "l",
                 xlab = "mcmc iterations",
                 ylab = "sampled values",
                 main = "complete trace (no burn-in)")
  graphics::abline(h = posterior_means[plot_num], col = "red")
  if (TRUE_VAL_AVAIL) {
    graphics::abline(h = true_vals[plot_num], col = "green")
  }
  graphics::mtext(par_names_plots[plot_num],
                  outer = TRUE,
                  cex = 1.5,
                  side = 3)
}
#' Generates ggplot graphic analysis plot of MCMC output.
#'
#' Generates ggplot graphic analysis plot of MCMC output for a model parameter
#' including: 1. histogram of mcmc draws (after burn-in and thinning) 2. trace
#' plot of mcmc draws (after burn-in and thinning) 3. autorcorrelation plot of
#' mcmc draws (after burn-in and thinning) 4. trace plot of mcmc draws (BEFORE
#' burn-in and thinning)
#'
#' @param mcmc_sims_df mcmc draws as a data.frame
#' @param mcmc_sims_df_after mcmc draws as a data.frame after burn-in and
#'   thinning
#' @param burn burn-in period
#' @param thin thinning; \code{thin = 10} means every 10th draw is excluded
#' @param num_mcmc total number of MCMC draws
#' @param par_names parameter names as used in the pgas output (col names for
#'   \code{mcmcm_sims} e.g.)
#' @param true_vals true values if a simulation is run; default is NULL and true
#'   values will not be added for the histogram and trace plots
#' @param posterior_means posterior means (have to be pre-computed and passed
#'   directly here)
#' @param plot_num number of parameter or plot
#'
#' @return ggplot in 2x2 form displaying the MCMC diagnostics
generate_ggplot2 <- function(mcmc_sims_df,
                             mcmc_sims_df_after,
                             burn = 0,
                             thin = 0,
                             num_mcmc,
                             par_names,
                             true_vals = NULL,
                             posterior_means,
                             plot_num) {
  TRUE_VAL_AVAIL <- !is.null(true_vals) && !is.na(true_vals)
  par_to_plot <- parse(text = par_names[plot_num])
  acfs    <- stats::acf(mcmc_sims_df_after[par_names[plot_num]], plot = FALSE)
  acfs_df <- with(acfs, data.frame(acfs$lag, acfs$acf))

  hist_plot <- ggplot2::ggplot(data = mcmc_sims_df_after) +
    ggplot2::geom_density(mapping = ggplot2::aes_string(x = par_names[plot_num])) +
    ggplot2::geom_histogram(ggplot2::aes(x = eval(par_to_plot),
                                         y = after_stat(density)),
                            binwidth = 0.025,
                            alpha = 0.5) +
    ggplot2::geom_vline(xintercept = posterior_means[plot_num], colour = "red") +
    ggplot2::xlab("values (after burn-in)") +
    ggplot2::ggtitle("posterior density histogram")
  if (TRUE_VAL_AVAIL) {
    hist_plot <- hist_plot +
      ggplot2::geom_vline(xintercept = true_vals[plot_num], colour = "green")
  }
  trace_plot_full <- ggplot2::ggplot(data = mcmc_sims_df,
                                     mapping = ggplot2::aes(x = num_mcmc,
                                                            y = eval(par_to_plot))) +
    ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept = posterior_means[plot_num], colour = "red") +
    ggplot2::ylab("sampled values") +
    ggplot2::xlab("mcmc iterations") +
    ggplot2::ggtitle("complete trace (no burn-in)")
  if (TRUE_VAL_AVAIL) {
    trace_plot_full <- trace_plot_full +
      ggplot2::geom_hline(yintercept = true_vals[plot_num], colour = "green")
  }

  trace_plot_burn <- ggplot2::ggplot(data = mcmc_sims_df_after,
                                     mapping = ggplot2::aes(x = num_mcmc,
                                                            y = eval(par_to_plot))) +
    ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept = posterior_means[plot_num], colour = "red") +
    ggplot2::ylab("sampled values") +
    ggplot2::xlab("mcmc iterations (after burn-in)") +
    ggplot2::ggtitle(paste("trace after burn-in", burn, sep = ": "))
  if (TRUE_VAL_AVAIL) {
    trace_plot_full <- trace_plot_full +
      ggplot2::geom_hline(yintercept = true_vals[plot_num], colour = "green")
  }

  acf_plot <- ggplot2::ggplot(data = acfs_df,
                              mapping = ggplot2::aes(x = .data$`acfs.lag`,
                                                     y = .data$`acfs.acf`)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
    ggplot2::geom_segment(mapping = ggplot2::aes(xend = .data$`acfs.lag`,
                                                 yend = 0)) +
    ggplot2::xlab("lag (after burn-in)") +
    ggplot2::ggtitle("autocorrelation")

  return(list(hist_plot, trace_plot_full, acf_plot, trace_plot_burn))
}
generate_ggplot2_all <- function(mcmc_sims_df,
                                 mcmc_sims_df_after,
                                 burn,
                                 thin,
                                 num_mcmc,
                                 par_names,
                                 true_vals,
                                 posterior_means,
                                 settings_plots,
                                 lab_names) {
  num_par <- length(par_names)
  list_of_plots <- vector("list", num_par)
  for (i in 1:num_par) {
    ln_tmp <- grid::textGrob(lab_names[i],
                             gp = grid::gpar(fontsize = 20, font = 3))
    plot_returned <- generate_ggplot2(mcmc_sims_df = mcmc_sims_df,
                                      mcmc_sims_df_after=mcmc_sims_df_after,
                                      burn = burn,
                                      thin = thin,
                                      num_mcmc = num_mcmc,
                                      par_names = par_names,
                                      true_vals = true_vals,
                                      posterior_means = posterior_means,
                                      plot_num = i)
    if (settings_plots$plot_view) {
      gridExtra::grid.arrange(plot_returned[[1]],
                              plot_returned[[2]],
                              plot_returned[[3]],
                              plot_returned[[4]],
                              nrow = 2,
                              top = ln_tmp)
    }
    if (settings_plots$plot_save) {
      current_plot_name <- get_plot_name(settings_plots,
                                         par_names[i],
                                         all = FALSE)
      list_of_plots[[i]] <- gridExtra::arrangeGrob(plot_returned[[1]],
                                                   plot_returned[[2]],
                                                   plot_returned[[3]],
                                                   plot_returned[[4]],
                                                   nrow = 2,
                                                   top = ln_tmp)
      ggplot2::ggsave(filename = basename(current_plot_name),
                      plot = list_of_plots[[i]],
                      path = dirname(current_plot_name),
                      device = "eps",
                      width = 18,
                      height = 10.5,
                      units = "cm")
      print(paste("Saved plots in: ", current_plot_name))
    }
    if (settings_plots$plot_save_all) {
      warning("Feature not yet implemented.")
    }
  }
}
generate_plot_all <- function(mcmc_sims,
                              mcmc_sims_after,
                              burn,
                              thin,
                              num_mcmc,
                              par_names,
                              true_vals,
                              posterior_means,
                              settings_plots,
                              lab_names) {
  num_par   <- length(par_names)
  if (settings_plots$plot_view) {
    for (i in 1:num_par) {
      generate_plot2(mcmc_sims = mcmc_sims,
                     mcmc_sims_after = mcmc_sims_after,
                     burn = burn,
                     num_mcmc = num_mcmc,
                     par_names = par_names,
                     par_names_plots = lab_names,
                     true_vals = true_vals,
                     posterior_means = posterior_means,
                     plot_num = i)

    }
  }
  if (settings_plots$plot_save) {
    for (i in 1:num_par) {
      current_plot_name <- get_plot_name(settings_plots,
                                         par_names[i],
                                         all = FALSE)
      grDevices::setEPS()
      grDevices::postscript(current_plot_name, width = 18,
                            height = 10.5)
      generate_plot2(mcmc_sims = mcmc_sims,
                     mcmc_sims_after = mcmc_sims_after,
                     burn = burn,
                     num_mcmc = num_mcmc,
                     par_names = par_names,
                     par_names_plots = lab_names,
                     true_vals = true_vals,
                     posterior_means = posterior_means,
                     plot_num = i)
      grDevices::dev.off()
      print(paste("Saved plots in: ", current_plot_name))
    }
  }
  if (settings_plots$plot_save_all) {
    current_plot_name_all <- get_plot_name(settings_plots, all = TRUE)
    grDevices::setEPS()
    grDevices::postscript(current_plot_name_all,
                          width = 18, height = 10.5,
                          onefile = TRUE)
    for (i in 1:num_par) {
      generate_plot2(mcmc_sims = mcmc_sims,
                     mcmc_sims_after = mcmc_sims_after,
                     burn = burn,
                     num_mcmc = num_mcmc,
                     par_names = par_names,
                     par_names_plots = lab_names,
                     true_vals = true_vals,
                     posterior_means = posterior_means,
                     plot_num = i)
      print(paste("Saving plots to big (.eps) file: ",
                  i, " out of ", num_par, "..."))
    }
    grDevices::dev.off()
  }
}
get_plot_name <- function(settings, sub_name, all) {
  if (all) {
    file.path(settings$plot_path,
              paste(settings$plot_name,
                    ".eps",
                    sep = ""))
  } else {
    file.path(settings$plot_path,
              paste(sub_name,
                    ".eps",
                    sep = ""))
  }

}
