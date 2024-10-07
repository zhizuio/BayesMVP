#' @title coef method for class \code{BayesMVP}
#' @description
#' Extract the posterior mean of the coefficients of a \code{BayesMVP} class object
#' @name coef.BayesMVP
#'
#' @param object an object of class \code{BayesMVP}
#' @param beta.type type of output beta. Default is \code{marginal}, giving
#' marginal beta estimation. If \code{beta.type="conditional"}, it gives beta
#' estimation conditional on gamma=1.
#' @param Pmax If \code{Pmax=0.5} and \code{beta.type="conditional"}, it gives
#' median probability model betas. Default is 0.
#' @param ... other arguments
#'
#' @return Estimated coefficients are from an object of class \code{BayesMVP}.
#' If the \code{BayesMVP} specified data standardization, the fitted values
#' are base based on standardized data.
#'
#' @examples
#' x <- 1
#'
#' @export
coef.BayesMVP <- function(object, beta.type = "marginal", Pmax = 0, ...) {
  if (!(beta.type %in% c("marginal", "conditional"))) {
    stop("Please specify correct 'beta.type'!")
  }
  if (Pmax < 0 || Pmax > 1) {
    stop("Please specify a correct argument 'Pmax' in [0,1]!")
  }
  if ((Pmax > 0) && (beta.type == "marginal")) {
    stop("Pmax > 0 is valid only if the argument beta.type = 'conditional'!")
  }
  getEstimator(object,
    estimator = "beta", Pmax = Pmax,
    beta.type = beta.type, ...
  )
}
