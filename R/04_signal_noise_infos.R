#' Information for diagnostics on regressors, parameters and measurements
#'
#' Under (P)MCMC sampling, convergence may deteriorate if regressors are higly
#' colinear, parameters are a-posteriori highly correlated or if the signal to
#' noise ratio is too high/low. This function helps to diagnose such problems by
#' printing:
#' \itemize{
#'   \item{regressor means: }{colmeans of regressor matrix}
#'   \item{regressor variance: }{col-variance of regressor matrix}
#'   \item{beta values: }{true parameter values for the betas}
#'   \item{mean of reg % * % beta: }{mean of this vector (matrix product)}
#'   \item{variance of reg % * % beta: }{variance of this vector (matrix product)}
#'   \item{mean of measurments Y: }{-}
#'   \item{variance of measurments Y: }{-}
#' }
#'
#'
#' @param regs1 first regressor matrix
#' @param regs2 second regressor matrix
#' @param beta1 first parameter vector
#' @param beta2 second parameter vector
#' @param Y1 first measurement component
#' @param Y2 second measurement component
#' @param names_compare a character string of two: names for comparisons
#' @param WIDE logical; if \code{TRUE}, returns a \code{2 x 7} summary matrix,
#'   else a \code{7 x 2} summary matrix
#' @param ROUND a numeric integer giving the digits to round the matrix outputs
#'
#' @return a summary matrix of measures as outlined in the "Description"
#' @export
get_signal_to_noise_infos <- function(regs1, regs2, beta1, beta2, Y1, Y2,
                                      names_compare = NULL, WIDE = FALSE,
                                      ROUND = 6) {
  if (is.null(names_compare)) names_compare <- c("first reg", "second reg")
  if (length(names_compare) != 2) stop("Arg 'names_compare' must have length 2")

  dim_tmp <- length(beta1)
  colmeans1 <- colMeans(regs1)
  colmeans2 <- colMeans(regs2)

  colVars1  <- apply(regs1, 2, var)
  colVars2  <- apply(regs2, 2, var)

  meanY1 <- mean(Y1)
  meanY2 <- mean(Y2)

  varY1 <- var(Y1)
  varY2 <- var(Y2)

  regBeta1 <- regs1 %*% beta1
  regBeta2 <- regs2 %*% beta2

  MeanRegBeta1 <- mean(regBeta1)
  MeanRegBeta2 <- mean(regBeta2)

  VarRegBeta1 <- var(regBeta1)
  VarRegBeta2 <- var(regBeta2)

  out <- matrix(c(colmeans1, colmeans2,
                  colVars1, colVars2,
                  beta1, beta2),
                ncol = 2*dim_tmp,
                byrow = TRUE)
  out <- round(out, digits = ROUND)
  out2 <- matrix(c(kappa(regs1), kappa(regs2),
                   kappa(t(regs1) %*% regs1),
                   kappa(t(regs2) %*% regs2)),
                 byrow = TRUE, ncol = 2)
  out3 <- matrix(c(MeanRegBeta1, VarRegBeta1, meanY1, varY1,
                   MeanRegBeta2, VarRegBeta2,
                   meanY2, varY2), ncol = 2)
  out3 <- round(out3, digits = ROUND)
  colnames(out)  <- rep(names_compare, each = dim_tmp)
  colnames(out2) <- names_compare
  colnames(out3) <- names_compare
  rownames(out)  <- c("regressor means: ",
                      "regressor variance: ",
                      "beta values: ")
  rownames(out2) <- c("kappa(regs)", "kappa(t(regs) %*% regs): ")
  rownames(out3) <- c("reg %*% beta means: ",
                      "reg %*% beta variance: ",
                      "Y means: ",
                      "Y variance: ")
  if (WIDE) return(list(t(out), t(out2), t(out3)))
  return(list(info1 = out,
              info2 = out2,
              info3 = out3))
}
