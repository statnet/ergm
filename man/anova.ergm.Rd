%  File ergm/man/anova.ergm.Rd
%  Part of the statnet package, http://statnet.org
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
%  Copyright 2012 the statnet development team
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{anova.ergm}
\title{ANOVA for Linear Model Fits}
\usage{
\method{anova}{ergm}(object, \dots, eval.loglik = FALSE)
\method{anova}{ergmlist}(object, \dots, eval.loglik = FALSE, scale = 0, test = "F")
}
\alias{anova.ergm}
\alias{anova.ergmlist}
\description{
  Compute an analysis of variance table for one or more linear model fits.
}
\arguments{
  \item{object, \dots}{objects of class \code{\link{ergm}}, usually, a result of a
    call to \code{\link{ergm}}.}
  \item{eval.loglik}{a logical specifying whether the log-likelihood will be
    evaluated if missing.}
  \item{test}{a character string specifying the test statistic to be
    used. Can be one of \code{"F"}, \code{"Chisq"} or \code{"Cp"},
    with partial matching allowed, or \code{NULL} for no test.}
  \item{scale}{numeric. An estimate of the noise variance
    \eqn{\sigma^2}{sigma^2}. If zero this will be estimated from the
    largest model considered.
  }
}
\details{
  Specifying a single object gives a sequential analysis of variance
  table for that fit.  That is, the reductions in the residual sum of
  squares as each term of the formula is added in turn are given in
  the rows of a table, plus the residual sum of squares.

  The table will contain F statistics (and P values) comparing the
  mean square for the row to the residual mean square.

  If more than one object is specified, the table has a row for the
  residual degrees of freedom and sum of squares for each model.
  For all but the first model, the change in degrees of freedom and sum
  of squares is also given. (This only make statistical sense if the
  models are nested.)  It is conventional to list the models from
  smallest to largest, but this is up to the user.

  Optionally the table can include test statistics.  Normally the
  F statistic is most appropriate, which compares the mean square for a
  row to the residual sum of squares for the largest model considered.
  If \code{scale} is specified chi-squared tests can be used. Mallows'
  \eqn{C_p}{Cp} statistic is the residual sum of squares plus twice the
  estimate of \eqn{\sigma^2}{sigma^2} times the residual degrees of
  freedom.

  If any of the objects do not have estimated log-likelihoods, produces
  an error, unless \code{eval.loglik=TRUE}.
}
\value{
  An object of class \code{"anova"} inheriting from class \code{"data.frame"}.
}
\section{Warning}{
  The comparison between two or more models will only be valid if they
  are fitted to the same dataset. This may be a problem if there are
  missing values and \R's default of \code{na.action = na.omit} is used,
  and \code{\link{anova.ergmlist}} will detect this with an error.
}
\seealso{
  The model fitting function \code{\link{ergm}}, \code{\link{anova}},
  \code{\link{logLik.ergm}} for adding the log-likelihood to an existing
  \code{\link[=ergm.object]{ergm}} object.
}
\examples{
data(molecule)
molecule \%v\% "atomic type" <- c(1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,3)
fit0 <- ergm(molecule ~ edges)
anova(fit0)
fit1 <- ergm(molecule ~ edges + nodefactor("atomic type"))
anova(fit1)
\donttest{
fit2 <- ergm(molecule ~ edges + nodefactor("atomic type") +  gwesp(0.5,
  fixed=TRUE), eval.loglik=TRUE) # Note the eval.loglik argument.
anova(fit0, fit1)
anova(fit0, fit1, fit2)
}
}
\keyword{regression}
\keyword{models}
