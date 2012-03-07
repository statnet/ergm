%  File ergm/man/coef.stergm.Rd
%  Part of the statnet package, http://statnet.org
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
%  Copyright 2012 the statnet development team
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{coef.stergm}
\alias{coef.stergm}
\alias{coefficients.stergm}
\title{Extract Model Coefficients}
\usage{
\method{coef}{stergm}(object, \dots)
\method{coefficients}{stergm}(object, \dots)
}
\arguments{
  \item{object}{A \code{\link{stergm}} fit.}
  \item{\dots}{other arguments.}
}
\value{
  Coefficients extracted from \code{object} in the form of a list with
  two elements: \verb{formation}, a vector of formation
  coefficients and \verb{dissolution} , a vector of dissolution
  coefficients.
}
\description{
  \code{coef} is a Method which extracts model coefficients
  from objects returned by the \code{\link{stergm}} function.  
  \code{coefficients} is an \emph{alias} for it.
}
\seealso{
  \code{\link{stergm}} and  \code{\link{ergm}}
}

\keyword{regression}
\keyword{models}
