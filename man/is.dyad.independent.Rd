%  File man/is.dyad.independent.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2014 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{is.dyad.independent}
\alias{is.dyad.independent}
\alias{is.dyad.independent.formula}
\alias{is.dyad.independent.ergm}
\title{
Testing for dyad-independence
}
\description{
These functions test whether an ERGM fit or formula is
dyad-independent.
}
\usage{
\method{is.dyad.independent}{ergm}(object, ...)
\method{is.dyad.independent}{formula}(object, 
                                      response=NULL,
                                      basis=NULL,
                                      \dots)
}
\arguments{
  \item{object}{
    An \code{\link{ergm}} object or an ERGM formula.
  }
  \item{response, basis}{
    See \code{\link{ergm}}.
  }
  \item{\dots}{
    Unused at this time.
  }
}
\value{
  \code{TRUE} if the model fit or one implied by the formula is
  dyad-independent; \code{FALSE} otherwise.
}

\keyword{ model }
