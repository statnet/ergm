%  File man/vcov.ergm.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2017 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{vcov.ergm}
\alias{vcov.ergm}
\title{Extract Model Covariance Matrix}
\usage{
\method{vcov}{ergm}(object, sources=c("all","model","estimation"), \dots)
}
\arguments{
  \item{object}{The \code{ergm} output object.}
  \item{sources}{Specify whether to return the covariance matrix from the ERGM model, 
  the estimation process, or both combined.}
  \item{\dots}{other arguments.}
}
\value{
  Coefficients extracted from the model object \code{object}.
}
\description{
  \code{vcov} is a method which extracts the covariance matrix
  from the output object returned by the \code{\link{ergm}} estimation.  
}
%\details{
%  All object classes which are returned by model fitting functions
%  should provide a \code{coef} method or use the default one.
%  (Note that the method is for \code{coef} and not \code{coefficients}.)
%
%  Class \code{"aov"} has a \code{coef} method that does not report
%  aliased coefficients (see \code{\link{alias}}).
%}
\seealso{
  \code{\link{coef.ergm}} 
}
%\references{
%  Chambers, J. M. and Hastie, T. J. (1992)
%  \emph{Statistical Models in S}.
%  Wadsworth \& Brooks/Cole.
%}
\examples{
data(florentine)
fit <- ergm(flomarriage ~ edges + concurrent)
vcov(fit, sources="model")
vcov(fit, sources="estimation")
vcov(fit, sources="all") # the default
}

