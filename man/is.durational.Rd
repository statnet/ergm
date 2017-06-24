%  File man/is.durational.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2017 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{is.durational}
\alias{is.durational}
\alias{is.durational.formula}
\alias{is.durational.ergm.model}
\alias{is.durational.NULL}
\alias{is.durational.character}
\title{
Testing for durational dependent models
}
\description{
These functions test whether an ERGM model or formula is
durational dependent or not. If the formula or model does not include any terms that need information about the duration of existing ties, the ergm proceass can use more efficient internal data structures.
}
\usage{
\method{is.durational}{character}(object, ...)
\method{is.durational}{ergm.model}(object, ...)
\method{is.durational}{formula}(object, 
                                      response=NULL,
                                      basis=NULL,
                                      \dots)
}
\arguments{
  \item{object}{
    An \code{\link{ergm}} object or an ERGM formula, or some characters, 
    e.g., object="all" for monitoring purpose.
  }
  \item{response, basis}{
    See \code{\link{ergm}}.
  }
  \item{\dots}{
    Unused at this time.
  }
}
\value{
  \code{TRUE} if the ERGM terms in the formula or model are durational dependent
  ; \code{FALSE} otherwise.
}

\keyword{ model }
