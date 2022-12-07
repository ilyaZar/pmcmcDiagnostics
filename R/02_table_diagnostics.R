#' Returns output diagnostics table
#'
#' The output diagnostics table summarizes posterior mean, standard deviation of
#' the parameter under the posterior distribution and standard deviation of the
#' posterior mean, confidence bands, HPDs, effective sample sizes (also as Stan
#' version with ESS bulk/tail) in columns for each parameter (in rows).
#'
#' @param num_par number of total parameters (rows of the output table)
#' @param par_names names of parameters to be taken as row-names in output table
#' @param mcmc_sims (particle) MCMC draws/simulated samples
#' @param mcmc_sims_after (particle) MCMC draws/simulated samples after burn-in
#' @param burn burn-in in period
#' @param num_mcmc number of MCMC iterations
#' @param posterior_means pre-computed posterior means
#' @param start_vals starting values of the parameter draws
#' @param true_vals true values with default set to \code{NULL}
#' @param ki_prob probability mass to cover with confidence bands and HPD
#' @param ESS_STANDARD logical; if \code{TRUE} computes the effective sample
#'   size
#' @param ESS_STAN logical; if \code{TRUE} computes the effective sample
#'   size in Stan-style i.e. with ESS bulk and ESS tail
#' @param settings_table settings for table output
#'
#' @return returns output diagnostic table with at least 8 columns (mean, sd,
#'   sd-posterior mean, confidence interval, HPDs) but up to 11 columns (if ESS,
#'   ESS bulk and ESS tail are included)
#' @export
diagnostics_table <- function(num_par,
                              par_names,
                              mcmc_sims,
                              mcmc_sims_after,
                              burn,
                              num_mcmc,
                              posterior_means,
                              start_vals,
                              true_vals = NULL,
                              ki_prob,
                              q_probs,
                              ESS_STANDARD = TRUE,
                              ESS_STAN = TRUE,
                              settings_table) {

  summary_results <- data.frame(start_val = numeric(num_par),
                                mean = numeric(num_par),
                                sd = numeric(num_par),
                                sd_mean = numeric(num_par),
                                CI_lo = numeric(num_par),
                                CI_up = numeric(num_par),
                                contained_CI = numeric(num_par),
                                HPD_lo = numeric(num_par),
                                HPD_up = numeric(num_par),
                                contained_HPD = numeric(num_par))
  if (!is.null(q_probs)) {
    q_probs_mat <- matrix(0, nrow = num_par, ncol = length(q_probs))
    summary_results <- cbind(summary_results, q_probs_mat)
    id_qs <- 11:(11 + length(q_probs) - 1)
    names(summary_results)[id_qs] <- paste0(q_probs * 100, "%")
  }
  if (ESS_STANDARD && !ESS_STAN) {
    summary_results = cbind(summary_results, ess = numeric(num_par))
  }
  if (ESS_STAN && !ESS_STANDARD) {
    summary_results = cbind(summary_results, ess_bulk = numeric(num_par))
    summary_results = cbind(summary_results, ess_tail = numeric(num_par))
    # The ess_bulk function produces an estimated Bulk Effective Sample Size
    # (bulk-ESS) using rank normalized draws. Bulk-ESS is useful measure for
    # sampling efficiency in the bulk of the distribution (related e.g. to
    # efficiency of mean and median estimates), and is well defined even if the
    # chains do not have finite mean or variance.
    #
    # The ess_tail function produces an estimated Tail Effective Sample Size
    # (tail-ESS) by computing the minimum of effective sample sizes for 5% and
    # 95% quantiles. Tail-ESS is useful measure for sampling efficiency in the
    # tails of the distribution (related e.g. to efficiency of variance and tail
    # quantile estimates).
    #
    # Both bulk-ESS and tail-ESS should be at least 100 (approximately) per
    # Markov Chain in order to be reliable and indicate that estimates of
    # respective posterior quantiles are reliable.
  }
  if (ESS_STANDARD && ESS_STAN) {
    summary_results = cbind(summary_results, ess = numeric(num_par))
    summary_results = cbind(summary_results, ess_bulk = numeric(num_par))
    summary_results = cbind(summary_results, ess_tail = numeric(num_par))
  }
  # summary_results2$true_vals <- true_vals
  cis <- compute_ci(mcmc_sims_after, ki_prob)
  hpd <- compute_hpd(mcmc_sims, mcmc_sims_after, burn, num_mcmc, ki_prob)
  ess <- compute_ess(mcmc_sims_after, ESS_STANDARD, ESS_STAN)
  qs  <- compute_quantiles(mcmc_sims_after, q_probs)
  summary_results$start_val <- start_vals
  summary_results$mean      <- posterior_means
  summary_results[, 3:4]    <- compute_sds(mcmc_sims_after, ess[, 1])
  summary_results[, 5:6]    <- cis
  summary_results[, 7]      <- compute_significance_indicator(cis, true_vals)
  summary_results[, 8:9]    <- hpd
  summary_results[, 10]     <- compute_significance_indicator(hpd, true_vals)
  summary_results[, id_qs]  <- qs
  id_ess <- (id_qs[length(id_qs)] + 1):ncol(summary_results)
  summary_results[, id_ess] <- compute_ess(mcmc_sims_after,
                                           ESS_STANDARD,
                                           ESS_STAN)
  if (!is.null(true_vals)) summary_results <- cbind(true_vals = true_vals,
                                                    summary_results)
  row.names(summary_results) <- unname(par_names)

  if (settings_table$table_view) {
    view_diagnostics_table(summary_results,
                           precision = settings_table$table_prec,
                           table_name = settings_table$table_name)
  }
  if (settings_table$table_save) {
    if (!is.null(settings_table$table_prec)) {
      summary_results <- round(summary_results,
                               digits = settings_table$table_prec)
    }
    write_diagnostics_table(summary_results, par_names,
                            settings_table$table_path,
                            paste0("SUMMARY_", settings_table$table_name))
  }
  ids <- compute_id_parts(summary_results, true_vals, q_probs)
  part_1 <- summary_results[, ids[["id_part1"]]]
  part_2 <- summary_results[, ids[["id_part2"]]]
  part_3 <- summary_results[, ids[["id_part3"]]]
  part_4 <- summary_results[, ids[["id_part4"]]]

  if(ncol(part_4) == 0) part_4 <- NULL

  return(list(all = summary_results,
              core = part_1,
              CIs = part_2,
              quantiles = part_3,
              convergence = part_4))
}
compute_id_parts <- function(summary_results, true_vals, q_probs) {
  ncol_sdmean <- which(names(summary_results) == "sd_mean")
  ncol_cnthpd <- which(names(summary_results) == "contained_HPD")
  ncol_qprobs <- ncol_cnthpd + length(q_probs)
  id_part1 <- seq_len(ncol_sdmean)
  id_part2 <- (ncol_sdmean + 1):ncol_cnthpd
  id_part3 <- (ncol_cnthpd + 1):ncol_qprobs

  if (ncol_qprobs < ncol(summary_results)) {
    id_part4 <- (ncol_qprobs + 1):ncol(summary_results)
  } else {
    id_part4 <- numeric(0)
  }
  return(list(id_part1 = id_part1, id_part2 = id_part2,
              id_part3 = id_part3, id_part4 = id_part4))
}
#' Compute HPD interval
#'
#' Computes before and after burn-in, and then chooses the smallest of those.
#'
#' @inheritParams diagnostics_table
#'
#' @return a two-column matrix that contains the smallest HPD interval for each
#'   parameter
#' @export
compute_hpd <- function(mcmc_sims, mcmc_sims_after,
                        burn, num_mcmc, ki_prob) {
  mcmc_sims_coda <- coda::mcmc(mcmc_sims, start = burn, end = num_mcmc)
  hpd_interval1  <- unname(coda::HPDinterval(mcmc_sims_coda, prob = ki_prob))
  hpd_interval2  <- t(HDInterval::hdi(mcmc_sims_after,
                                      credMass = ki_prob,
                                      allowSplit = TRUE))
  min_hpd <- cbind(HPD1 = hpd_interval1[, 2] - hpd_interval1[, 1],
                   HPD2 = hpd_interval2[, 2] - hpd_interval2[, 1])
  min_hpd <- apply(min_hpd, 1, which.min)

  num_pars <-  nrow(hpd_interval1)
  hpd_interval_out <- matrix(0, nrow = num_pars,
                             ncol = ncol(hpd_interval1))
  for (i in 1:num_pars) {
    if (min_hpd[i] == 1) {
      hpd_interval_out[i, ] <- hpd_interval1[i, ]
    } else {
      hpd_interval_out[i, ] <- hpd_interval2[i, ]
    }
  }
  colnames(hpd_interval_out) <- c("HPD_lo", "HPD_up")
  return(hpd_interval_out)
}
#' Computes the expected sample size.
#'
#' Uses samples after burn-in period.
#'
#' @inheritParams diagnostics_table
#'
#' @return a vector of expected sample sizes of length equal to the number of
#'   parameters
#' @export
compute_ess <- function(mcmc_sims_after, ESS_STANDARD, ESS_STAN) {
  num_pars <- ncol(mcmc_sims_after)
  if (ESS_STANDARD && !ESS_STAN) {
    ess <- numeric(num_pars)
    for (i in 1:num_pars) {
      ess[i] <- tryCatch(mcmcse::ess(mcmc_sims_after[, i]),
                         error = function(err) NA_real_)
    }
    ess <- matrix(round(ess, digits = 0), ncol = 1)
    colnames(ess) <- "ess"
  } else if (ESS_STAN && !ESS_STANDARD) {
    ess_bulk <- numeric(num_pars)
    ess_tail <- numeric(num_pars)
    for (i in 1:num_pars) {
      ess_bulk[i] <- tryCatch(rstan::ess_bulk(mcmc_sims_after[, i]),
                              error = function(err) NA_real_)
      ess_tail[i] <- tryCatch(rstan::ess_tail(mcmc_sims_after[, i]),
                              error = function(err) NA_real_)
    }
    ess_bulk <- round(ess_bulk, digits = 0)
    ess_tail <- round(ess_bulk, digits = 0)
    ess <- cbind(ess_bulk = ess_bulk, ess_tail = ess_tail)
  } else if (ESS_STAN && ESS_STANDARD) {
    ess_stnd <- numeric(num_pars)
    ess_bulk <- numeric(num_pars)
    ess_tail <- numeric(num_pars)
    for (i in 1:num_pars) {
      ess_bulk[i] <- tryCatch(rstan::ess_bulk(mcmc_sims_after[, i]),
                              error = function(err) NA_real_)
      ess_tail[i] <- tryCatch(rstan::ess_tail(mcmc_sims_after[, i]),
                              error = function(err) NA_real_)
      ess_stnd[i] <- tryCatch(mcmcse::ess(mcmc_sims_after[, i]),
                              error = function(err) NA_real_)
    }
    ess_stnd <- round(ess_stnd, digits = 0)
    ess_bulk <- round(ess_bulk, digits = 0)
    ess_tail <- round(ess_tail, digits = 0)
    ess <- cbind(ess_stnd = ess_stnd,
                 ess_bulk = ess_bulk,
                 ess_tail = ess_tail)
  }
  return(ess)
}
#' Computes posterior standard deviation & standard deviation of posterior mean
#'
#' Standard deviation of the posterior mean is based on the expected sample
#' size.
#'
#' @inheritParams diagnostics_table
#' @param ess expected sample sizes as produced by [compute_ess]
#'
#' @return a two column matrix of two types of standard deviations:
#'   \itemize{
#'     \item posterior standard deviation
#'     \item estimate of standard deviation of the posterior mean based on the
#'     expected sample size
#'   }
#' @export
compute_sds <- function(mcmc_sims_after, ess) {
  num_pars <- length(ess)

  out1 <- numeric(num_pars)
  out2 <- numeric(num_pars)
  for(i in 1:num_pars) {
    out1[i] <- stats::sd(mcmc_sims_after[, i])
    out2[i] <- out1[i]/sqrt(ess[i])
  }

  out <- cbind(out1, out2)
  colnames(out) <- c("sd", "sd_mean")
  return(out)
}
#' Computes confidence interval
#'
#' Confidence interval is based on samples after burn-in.
#'
#' @inheritParams diagnostics_table
#'
#' @return a two column matrix with upper/lower (first/second column) bound of
#'   a confidence interval with significance given by \code{ki_prob}
#' @export
compute_ci <- function(mcmc_sims_after, ki_prob) {
  num_pars <- ncol(mcmc_sims_after)

  CI <- matrix(0, nrow = num_pars, ncol = 2)
  for (i in 1:num_pars) {
    CI[i, ] <- stats::quantile(mcmc_sims_after[, i],
                               probs = c((1 - ki_prob)/2, 1 - (1 - ki_prob)/2),
                               names = FALSE)
  }
  colnames(CI) <- c("CI_lo", "CI_up")
  return(CI)
}
compute_quantiles <- function(mcmc_sims_after,
                              q_probs) {
  num_pars <- ncol(mcmc_sims_after)

  quants <- matrix(0, nrow = num_pars, ncol = length(q_probs))
  for (i in 1:num_pars) {
    quants[i, ] <- stats::quantile(mcmc_sims_after[, i],
                                   probs = q_probs,
                                   names = FALSE)
  }
  colnames(quants) <- paste0(q_probs * 100, "%")
  return(quants)
}
#' Checks whether HPD is smaller than CI (for the same significance)
#'
#' Ideally, the HPD is the smallest/densest region around the posterior mean,
#' so there are computational issues whenever this is not the case (and, here
#' specifically, whenever the HPD is not smaller than the confidence interval
#' at the same significance level).
#'
#' @param CI a confidence interval as returned via [compute_ci]
#' @param HPD a confidence interval as returned via [compute_hpd]
#'
#' @return side effect: message indicating when HPD is not smaller than CI
#' @export
verify_CIs <- function(CI, HPD) {
  verify_CIs         <- cbind(CI[, 2] - CI[, 1], HPD[, 2] - HPD[, 1])
  check_hpd          <- apply(verify_CIs, 1, which.min)
  id_check_phd_fails <- which(check_hpd == 1)
  check_hpd          <- unique(check_hpd)
  if (!(length(check_hpd) == 1) || !(check_hpd == 2)) {
    msg1 <- paste0("HPD is not smaller than CI intervall for param. numbers: ")
    msg2 <- paste0(as.character(id_check_phd_fails), collapse = ", ")
    msg3 <- "\n"
    cat(crayon::red(paste0(msg1, msg2, msg3)))
  }
  return(invisible(NULL))
}
#' Check whether CI and HPD contain true values or zeros.
#'
#' The former case is for simulations (with \code{true_vals} being not
#' \code{NULL}) and the latter checks for significance i.e. if 0 is covered or
#' not.
#'
#' @param int a matrix of two columns with upper/lower confidence bands as
#'   produced by [compute_ci] or [compute_hpd]
#' @inheritParams diagnostics_table
#'
#' @return a list of two vectors giving \code{TRUE} if \code{true_vals} is
#'   covered (or zero, whenever true values are not passed), or \code{FALSE},
#'   if that is not the case for both CI and HPD
#' @export
compute_significance_indicator <- function(int, true_vals) {
  if (!is.null(true_vals)) {
    contained_CI  <-  (int[, 1] <= true_vals) & (int[, 2] >= true_vals)
  } else {
    contained_CI  <- sign(int[, 1]) == sign(int[, 2])
  }
  return(contained_CI)
}
#' Call \code{View} on \code{data.frame} of MCMC convergence diagnostics
#'
#' @param summary_diagnostics \code{data.frame} of MCMC convergence diagnostics
#' @param precision digits used for rounding values
#' @param table_name name of the table
#'
#' @return pure side effect call to \code{View}
#' @export
view_diagnostics_table <- function(summary_diagnostics, precision, table_name) {
  ID_round <- which(!sapply(summary_diagnostics, is.logical))
  summary_diagnostics[, ID_round] <- round(summary_diagnostics[, ID_round],
                                           digits = precision)
  utils::View(summary_diagnostics,
              title = paste(table_name, "_summary_diagnostics", sep = ""))
}
#' Write \code{data.frame} of convergence diagnostics to '.csv'
#'
#' @param summary_diagnostics \code{data.frame} of MCMC convergence diagnostics
#' @param par_names parameter names
#' @param table_path path to write to
#' @param table_name name of the table
#'
#' @return pure side effect; writing summary table to \code{.csv} file
#' @export
write_diagnostics_table <- function(summary_diagnostics, par_names,
                                    table_path, table_name) {
  if (is.null(par_names)) {
    stop("Can't save results in table form: label names required!")
  }
  summary_diagnostics_save <- summary_diagnostics
  summary_diagnostics_save <- cbind(par_names, summary_diagnostics_save)
  # row.names(summary_diagnostics_save) <- par_names
  # summary_diagnostics_save <- cbind(par_name = par_names, summary_diagnostics_save)
  write.csv(summary_diagnostics_save,
            file = file.path(table_path,paste0(table_name, ".csv")),
            row.names = FALSE)
}
