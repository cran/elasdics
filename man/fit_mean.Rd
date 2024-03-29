% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_mean.R
\name{fit_mean}
\alias{fit_mean}
\title{Fitting function for open curves}
\usage{
fit_mean(srv_data_curves, knots, max_iter, type, eps)
}
\arguments{
\item{srv_data_curves}{list of \code{data.frame}s with srv vectors in each row.
Usually a result of a call to \code{\link{get_srv_from_points}}}

\item{knots}{set of knots for the mean spline curve}

\item{max_iter}{maximal number of iterations}

\item{type}{if "smooth" linear srv-splines are used which results in a differentiable mean curve
if "polygon" the mean will be piecewise linear.}

\item{eps}{the algorithm stops if L2 norm of coefficients changes less}
}
\value{
a \code{list}
with entries
  \item{type}{"smooth" or "polygon"}
  \item{coefs}{\code{coefs} srv spline coefficients of the estimated mean}
  \item{knots}{spline knots}
  \item{t_optims}{optimal parametrization}
}
\description{
Fits an elastic mean for open curves. Is usually called from
\code{\link{compute_elastic_mean}}.
}
