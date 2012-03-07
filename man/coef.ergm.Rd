%  File ergm/man/coef.ergm.Rd
%  Part of the statnet package, http://statnet.org
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
%  Copyright 2012 the statnet development team
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{coef.ergm}
\alias{coef.ergm}
\alias{coefficients.ergm}
\title{Extract Model Coefficients}
\usage{
\method{coef}{ergm}(object, \dots)
\method{coefficients}{ergm}(object, \dots)
}
\alias{coef.ergm}
\alias{coefficients.ergm}
\arguments{
  \item{object}{an object for which the extraction of model coefficients is
    meaningful.}
  \item{\dots}{other arguments.}
}
\value{
  Coefficients extracted from the model object \code{object}.
}
\description{
  \code{coef} is a Method which extracts model coefficients
  from objects returned by the \code{\link{ergm}} function.  
  \code{coefficients} is an \emph{alias} for it.
}
\seealso{
  \code{\link{fitted.values}} and \code{\link{residuals}} for related methods;
  \code{\link{glm}}, \code{\link{lm}} for model fitting.
}
\examples{
data(molecule)
molecule \%v\% "atomic type" <- c(1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,3)
fit <- ergm(molecule ~ edges + nodefactor("atomic type"))
coef(fit)
}
\keyword{regression}
\keyword{models}
