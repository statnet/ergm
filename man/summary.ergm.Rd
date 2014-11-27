%  File man/summary.ergm.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2014 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{summary.ergm}
\alias{summary.ergm}
\alias{print.summary.ergm}
%\alias{print.summary.lm}
% \alias{print.summary.lm.null}
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
% and additionally gives
%\dQuote{significance stars}. if \code{signif.stars} is \code{TRUE}.

% Correlations are printed to two decimal places (or symbolically): to
% see the actual correlations print \code{summary(object)$correlation}
% directly.
}
\value{
  The function \code{\link{summary.ergm}} computes and returns a list of summary
  statistics of the fitted \code{\link{ergm}} model given in \code{object}.
%, using
%  the components (list elements) \code{"call"} and \code{"terms"}
%  from its argument, plus
%  \item{residuals}{the \emph{weighted} residuals, the usual residuals
%    rescaled by the square root of the weights specified in the call to
%    \code{\link{ergm}}.}
%  \item{coefficients}{a \eqn{p \times 4}{p x 4} matrix with columns for
%    the estimated coefficient, its standard error, t-statistic and
%    corresponding (two-sided) p-value.  Aliased coefficients are omitted.}
%  \item{aliased}{named logical vector showing if the original
%    coefficients are aliased.} 
%  \item{sigma}{the square root of the estimated variance of the random
%    error
%    \deqn{\hat\sigma^2 = \frac{1}{n-p}\sum_i{R_i^2},}{%
%      sigma^2 = 1/(n-p) Sum(R[i]^2),}
%    where \eqn{R_i}{R[i]} is the \eqn{i}-th residual, \code{residuals[i]}.}
%  \item{df}{degrees of freedom, a 3-vector \eqn{(p, n-p, p*)}, the last
%    being the number of non-aliased coefficients.}
%  \item{fstatistic}{(for models including non-intercept terms)
%    a 3-vector with the value of the F-statistic with
%    its numerator and denominator degrees of freedom.}
%  \item{r.squared}{\eqn{R^2}, the \dQuote{fraction of variance explained by
%    the model},
%    \deqn{R^2 = 1 - \frac{\sum_i{R_i^2}}{\sum_i(y_i- y^*)^2},}{%
%      R^2 = 1 - Sum(R[i]^2) / Sum((y[i]- y*)^2),}
%    where \eqn{y^*}{y*} is the mean of \eqn{y_i}{y[i]} if there is an
%    intercept and zero otherwise.}
%  \item{adj.r.squared}{the above \eqn{R^2} statistic
%    \dQuote{\emph{adjusted}}, penalizing for higher \eqn{p}.}
%  \item{cov.unscaled}{a \eqn{p \times p}{p x p} matrix of (unscaled)
%    covariances of the \eqn{\hat\beta_j}{coef[j]}, \eqn{j=1, \dots, p}.}
%  \item{correlation}{the correlation matrix corresponding to the above
%    \code{cov.unscaled}, if \code{correlation = TRUE} is specified.}
%  \item{symbolic.cor}{(only if \code{correlation} is true.)  The value
%    of the argument \code{symbolic.cor}.}
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
