%  File ergm/man/print.ergm.Rd
%  Part of the statnet package, http://statnet.org
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
%  Copyright 2012 the statnet development team
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{print.ergm}
\alias{print.ergm}
\title{Exponential Random Graph Models}
\description{
  \code{\link{print.ergm}} is the method used to print an \code{\link{ergm}} object created 
  by the \code{\link{ergm}} function.
}
\usage{
\method{print}{ergm}(x, digits = max(3, getOption("digits") - 3), \dots)
}
\arguments{
  \item{x}{An \code{\link{ergm}} object.  See documentation for \code{\link{ergm}}.}
  \item{digits}{Significant digits for coefficients}
  \item{\dots}{Additional arguments, to be passed to lower-level functions
  in the future.}
}
\details{Automatically called when an object of class \code{\link{ergm}} is printed.
  Currently, \code{\link{print.ergm}} summarizes 
  the size of the MCMC sample, the theta vector
  governing the selection of the sample, and the Monte Carlo MLE.
}
\value{The value returned is the \code{\link{ergm}} object itself.}
\seealso{network, ergm}
\examples{
 data(florentine)

 x <- ergm(flomarriage ~ density)
 class(x)
 x
}
\keyword{models}



