%  File ergm/man/print.stergm.Rd
%  Part of the statnet package, http://statnet.org
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
%  Copyright 2012 the statnet development team
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{print.stergm}
\alias{print.stergm}
\alias{print.summary.stergm}
\alias{summary.stergm}
\title{Separable Temporal Exponential Random Graph Models}
\description{
  \code{\link{print.stergm}} and \code{\link{summary.stergm}} are the method 
    used to print and summarize objects created by the \code{\link{stergm}}
    function.
}
\usage{
\method{print}{stergm}(x, digits = max(3, getOption("digits") - 3), \dots)
\method{summary}{stergm}(object, \dots)
}
\arguments{
  \item{x, object}{A \code{\link{stergm}} object.  See documentation for \code{\link{stergm}}.}
  \item{digits}{Significant digits for coefficients}
  \item{\dots}{Additional arguments, to be passed to lower-level functions
  in the future.}
}
\details{Automatically called when an object of class \code{\link{stergm}} is printed.
}
\value{The value returned is the \code{\link{stergm}} object itself.}
\seealso{\code{\link{network}}, \code{\link{ergm}}}
\keyword{models}



