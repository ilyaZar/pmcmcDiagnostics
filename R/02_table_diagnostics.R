#' Returns output diagnostics table
#'
#' The output diagnostics table summarizes posterior mean, standard deviation of
#' the parameter under the posterior distribution and standard deviation of the
#' posterior mean, confidence bands, HPDs, effective sample sizes (also as stan
#' version with ESS bulk/tail) in columns for each parameter (in rows).
#'
#' @param num_par number of total parameters (rows of the output table)
#' @param mcmc_sims (particle) MCMC draws/simulated samples
#' @param mcmc_sims_after (particle) MCMC draws/simulated samples after burnin
#' @param burn burnin in period
#' @param num_mcmc number of MCMC iterations
#' @param posterior_means precomputed posterior means
#' @param start_vals starting values of the parameter draws
#' @param true_vals true values with default set to \code{NULL}
#' @param ki_prob probability mass to cover with confidence bands and HPD
#' @param compute_ess logical; if \code{TRUE} computes the effective sample size
#' @param compute_ess_stan logical; if \code{TRUE} computes the effective sample
#'   size in stan-style i.e. with ESS bulk and ESS tail
#'
#' @return returns output diagnostic table with at least 8 columns (mean, sd,
#'   sd-posterior mean, confidence intervall, HPDs) but up to 11 columns (if
#'   ESS, ESS bulk and ESS tail are included)
#' @export
diagnostics_table <- function(num_par,
                              mcmc_sims,
                              mcmc_sims_after,
                              burn,
                              num_mcmc,
                              posterior_means,
                              start_vals,
                              true_vals = NULL,
                              ki_prob,
                              compute_ess = TRUE,
                              compute_ess_stan = TRUE) {
  summary_results <- data.frame(start_val = numeric(num_par),
                                mean = numeric(num_par),
                                sd = numeric(num_par),
                                sd_mean = numeric(num_par),
                                KI_lo = numeric(num_par),
                                KI_up = numeric(num_par),
                                HPD_lo = numeric(num_par),
                                HPD_up = numeric(num_par))
  if (compute_ess) {
    summary_results = cbind(summary_results, ess = numeric(num_par))
  }
  if (compute_ess_stan) {
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
  mcmc_sims_coda <- coda::mcmc(mcmc_sims, start = burn, end = num_mcmc)
  hpd_interval1  <- unname(coda::HPDinterval(mcmc_sims_coda, prob = ki_prob))
  hpd_interval2  <- t(HDInterval::hdi(mcmc_sims_after,
                                      credMass = ki_prob,
                                      allowSplit = TRUE))
  min_hpd <- cbind(HPD1 = hpd_interval1[, 2] - hpd_interval1[, 1],
                   HPD2 = hpd_interval2[, 2] - hpd_interval2[, 1])
  min_hpd <- apply(min_hpd, 1, which.min)
  ess <- numeric(num_par)
  for (i in 1:num_par) {
    val <- tryCatch(mcmcse::ess(mcmc_sims_after[, i]),
                    error = function(err) NA_real_)
    ess[i] <- val
  }
  for (i in 1:num_par) {
    summary_results[i, 1] <- start_vals[i]
    summary_results[i, 2] <- posterior_means[i]
    summary_results[i, 3] <- stats::sd(mcmc_sims_after[, i])
    summary_results[i, 4] <- summary_results[i, 3]/sqrt(ess[i])
    KI <- stats::quantile(mcmc_sims_after[, i],
                          probs = c((1 - ki_prob)/2, 1 - (1 - ki_prob)/2),
                          names = FALSE)
    summary_results[i, 5] <- KI[1]
    summary_results[i, 6] <- KI[2]
    if (min_hpd[i] == 1) {
      summary_results[i, 7] <- hpd_interval1[i, 1]
      summary_results[i, 8] <- hpd_interval1[i, 2]
    } else {
      summary_results[i, 7] <- hpd_interval2[i, 1]
      summary_results[i, 8] <- hpd_interval2[i, 2]
    }
    if (compute_ess) {
      summary_results[i, 9] <- round(ess[i], digits = 0)
    }
    if (compute_ess_stan && compute_ess) {
      summary_results[i, 10] <- round(rstan::ess_bulk(mcmc_sims_after[, i]),
                                      digits = 0)
      summary_results[i, 11] <- round(rstan::ess_tail(mcmc_sims_after[, i]),
                                      digits = 0)
    } else if (compute_ess_stan && !compute_ess) {
      summary_results[i, 9] <- round(rstan::ess_bulk(mcmc_sims_after[, i]),
                                     digits = 0)
      summary_results[i, 10] <- round(rstan::ess_tail(mcmc_sims_after[, i]),
                                      digits = 0)
    }
  }
  verify_KIs <- cbind(KI  = summary_results[, 6] - summary_results[, 5],
                      HPD = summary_results[, 8] - summary_results[, 7])
  check_hpd          <- apply(verify_KIs, 1, which.min)
  id_check_phd_fails <- which(check_hpd == 1)
  check_hpd          <- unique(check_hpd)
  if (!(length(check_hpd) == 1) || !(check_hpd == 2)) {
    msg1 <- paste0("HPD is not smaller than KI intervall for param. numbers: ")
    msg2 <- paste0(as.character(id_check_phd_fails), collapse = ", ")
    msg3 <- "\n"
    cat(crayon::red(paste0(msg1, msg2, msg3)))
  }
  if (!is.null(true_vals)) {
    contained_KI <-  ((summary_results[, 5, drop = TRUE] <= true_vals) &
                        (summary_results[, 6, drop = TRUE] >= true_vals))
    contained_hpd <-  ((summary_results[, 7, drop = TRUE] <= true_vals) &
                         (summary_results[, 8, drop = TRUE] >= true_vals))
    summary_results_true <- cbind(true_vals = true_vals,
                                  summary_results[, 1:6], contained_KI,
                                  summary_results[, 7:8], contained_hpd)
    num_add_rest <- ncol(summary_results) + 3 - ncol(summary_results_true)
    if (num_add_rest != 0) {
      summary_results_true <- cbind(summary_results_true,
                                    summary_results[, (8 + 1):(8+num_add_rest),
                                                    drop = FALSE])
    }
    return(summary_results_true)
  } else {
    significant_KI <-  (sign(summary_results[, 5, drop = TRUE]) ==
                        sign(summary_results[, 6, drop = TRUE]))
    significant_hpd <-  (sign(summary_results[, 7, drop = TRUE]) ==
                         sign(summary_results[, 8, drop = TRUE]))
    summary_results_significant <- cbind(summary_results[, 1:6],
                                         significant_KI,
                                         summary_results[, 7:8],
                                         significant_hpd)
    num_add_rest <- ncol(summary_results) + 2 -ncol(summary_results_significant)
    if (num_add_rest != 0) {
      summary_results_significant <- cbind(summary_results_significant,
                                    summary_results[, (8 + 1):(8 +num_add_rest),
                                                    drop = FALSE])
    }
    return(summary_results_significant)
  }
}
