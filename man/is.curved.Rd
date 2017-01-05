%  File man/is.curved.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2017 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{is.curved}
\alias{is.curved}
\alias{is.curved.formula}
\alias{is.curved.ergm}
\alias{is.curved.NULL}
\alias{is.curved.ergm.model}
\title{
Testing for curved exponential family
}
\description{
These functions test whether an ERGM fit or formula is
curved.
}
\usage{
\method{is.curved}{ergm}(object, ...)
\method{is.curved}{formula}(object, 
                    response=NULL,
                    basis=NULL,
                    \dots)
\method{is.curved}{ergm.model}(object, ...)
\method{is.curved}{NULL}(object, ...)
}
\arguments{
  \item{object}{
    An \code{\link{ergm}} object or an ERGM formula.
  }
  \item{response, basis}{
    Name of the edge attribute whose value is to be modeled. See \code{\link{ergm}}.
  }
  \item{\dots}{
    Unused at this time.
  }
}
\details{
Curvature is checked by testing if all model parameters are canonical. 
}
\value{
  \code{TRUE} if the model fit or one implied by the formula is
  curved; \code{FALSE} otherwise.
}

\keyword{ model }
