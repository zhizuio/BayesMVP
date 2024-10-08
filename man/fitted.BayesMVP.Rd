% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitted.BayesMVP.R
\name{fitted.BayesMVP}
\alias{fitted.BayesMVP}
\title{get fitted responses}
\usage{
\method{fitted}{BayesMVP}(object, Pmax = 0, beta.type = "marginal", ...)
}
\arguments{
\item{object}{an object of class \code{BayesMVP}}

\item{Pmax}{valid if \code{beta.type="conditional"}. If
\code{beta.type="conditional"} and \code{Pmax=0.5}, it gives median
probability model betas. Default is 0}

\item{beta.type}{type of estimated beta for the fitted model. Default is
\code{marginal}, giving marginal beta estimation. If
\code{beta.type="conditional"}, it gives beta estimation conditional
on gamma=1}

\item{...}{other arguments}
}
\value{
Fitted values extracted from an object of class \code{BayesMVP}. If
the \code{BayesMVP} specified data standardization, the fitted values are
base based on standardized data.
}
\description{
Return the fitted response values that correspond to the posterior mean
estimates from a \code{BayesMVP} class object.
}
\examples{
x <- 1

}
