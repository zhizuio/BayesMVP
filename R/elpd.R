#' @title expected log pointwise predictive density
#' @description
#' Measure the prediction accuracy by the elpd (expected log pointwise
#' predictive density). The out-of-sample predictive fit can either be
#' estimated by Bayesian leave-one-out cross-validation (LOO) or by widely
#' applicable information criterion (WAIC) (Vehtari et al. 2017).
#' @name elpd
#'
#' @param object an object of class \code{BayesMVP}
#' @param method the name of the prediction accuracy index. Default is the
#' \code{"LOO"} (Bayesian LOO estimate of out-of-sample predictive fit). The
#' other index is the \code{"WAIC"} (widely applicable information criterion).
#' For the HRR models, both "\code{LOO}" and "\code{WAIC}" are computed based
#' on the multivate t-distribution of the posterior predictive rather than
#' approximation of importance sampling.
#'
#' @return Return the predictiion accuracy measure from an object of class
#' \code{BayesMVP}. It is elpd.loo if the argumnet \code{method="LOO"} and
#' elpd.WAIC if \code{method="WAIC"}.
#'
#' @references Vehtari, A., Gelman, A., Gabry, J. (2017). \emph{Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC.} Statistics and Computing, 27(5): 1413–1432.
#'
#' @examples
#' x <- 1
#'
#' @export
elpd <- function(object, method = "LOO") {
  object$output[-1] <-
    paste(object$output$outFilePath, object$output[-1], sep = "")

  if (is.null(object$output$CPO)) {
    stop("Please specify the argument 'output_CPO = TRUE' in BayesMVP()!")
  }

  if (toupper(method) == "LOO") {
    elpd <- sum(log(read.table(object$output$CPO)))
    names(elpd) <- "elpd.loo"
  } else if (toupper(method) == "WAIC") {
    elpd <- sum(read.table(object$output$WAIC))
    names(elpd) <- "elpd.waic"
  } else {
    stop("Please give the correct method name!")
  }

  return(elpd)
}
