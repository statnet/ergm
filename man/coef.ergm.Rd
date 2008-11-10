%  File ergm/man/coef.ergm.Rd
%  Part of the statnet package, http://statnetproject.org
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnetproject.org/attribution
%
% Copyright 2003 Mark S. Handcock, University of Washington
%                David R. Hunter, Penn State University
%                Carter T. Butts, University of California - Irvine
%                Steven M. Goodreau, University of Washington
%                Martina Morris, University of Washington
% Copyright 2007 The statnet Development Team
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
%\details{
%  All object classes which are returned by model fitting functions
%  should provide a \code{coef} method or use the default one.
%  (Note that the method is for \code{coef} and not \code{coefficients}.)
%
%  Class \code{"aov"} has a \code{coef} method that does not report
%  aliased coefficients (see \code{\link{alias}}).
%}
\seealso{
  \code{\link{fitted.values}} and \code{\link{residuals}} for related methods;
  \code{\link{glm}}, \code{\link{lm}} for model fitting.
}
%\references{
%  Chambers, J. M. and Hastie, T. J. (1992)
%  \emph{Statistical Models in S}.
%  Wadsworth \& Brooks/Cole.
%}
\examples{
data(molecule)
molecule \%v\% "atomic type" <- c(1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,3)
fit <- ergm(molecule ~ edges + nodefactor("atomic type"))
coef(fit)
}
\keyword{regression}
\keyword{models}
