%  File ergm/man/summary.ergm.Rd
%  Part of the statnet package, http://statnet.org
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
%  Copyright 2012 the statnet development team
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{summary.ergm}
\alias{summary.ergm}
\alias{print.summary.ergm}
\title{Summarizing ERGM Model Fits}
\usage{
\method{summary}{ergm}(object, \dots,
               digits = max(3, getOption("digits") - 3),
               correlation = FALSE, covariance = FALSE, 
               total.variation=TRUE,
               eps = 1e-04)

}
\arguments{
  \item{object}{an object of class \code{"ergm"}, usually, a result of a
    call to \code{\link{ergm}}.}
  \item{digits}{Significant digits for coefficients}
  \item{correlation}{logical; if \code{TRUE}, the correlation matrix of
    the estimated parameters is returned and printed.}
  \item{covariance}{logical; if \code{TRUE}, the covariance matrix of
    the estimated parameters is returned and printed.}
  \item{total.variation}{logical; if \code{TRUE}, the standard errors
reported in the \code{Std. Error} column are based on the sum of the likelihood variation and the MCMC variation. If \code{FALSE} only the likelihood varuation is used. The \eqn{p}-values are based on this source of variation.}
  \item{eps}{number; indicates the smallest p-value. See \code{\link{printCoefmat}}.}
  \item{\dots}{Arguments to \code{\link{logLik.ergm}}}
}
\description{
\code{\link[base]{summary}} method for class \code{"ergm"}.
}
\details{
 \code{\link{summary.ergm}} tries to be smart about formatting the
 coefficients, standard errors, etc.

}
\value{
  The function \code{\link{summary.ergm}} computes and returns a list of summary
  statistics of the fitted \code{\link{ergm}} model given in \code{object}.
}
\seealso{ network, ergm, print.ergm.
  The model fitting function \code{\link{ergm}}, \code{\link{summary}}.

  Function \code{\link{coef}} will extract the matrix of coefficients
  with standard errors, t-statistics and p-values.
}
\examples{
 data(florentine)

 x <- ergm(flomarriage ~ density)
 summary(x)
}
\keyword{regression}
\keyword{models}
