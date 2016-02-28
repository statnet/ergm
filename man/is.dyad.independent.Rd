%  File man/is.dyad.independent.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2015 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{is.dyad.independent}
\alias{is.dyad.independent}
\alias{is.dyad.independent.formula}
\alias{is.dyad.independent.ergm}
\alias{is.dyad.independent.NULL}
\alias{is.dyad.independent.conlist}
\alias{is.dyad.independent.ergm.model}
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
\method{is.dyad.independent}{conlist}(object, 
                    object.obs = NULL, 
                    ...)
\method{is.dyad.independent}{ergm.model}(object, ...)
\method{is.dyad.independent}{NULL}(object, ...)
}
\arguments{
  \item{object}{
    An \code{\link{ergm}} object or an ERGM formula.
  }
  \item{response, basis}{
    Name of the edge attribute whose value is to be modeled. See \code{\link{ergm}}.
  }
  \item{object.obs}{
  
  }
  \item{\dots}{
    Unused at this time.
  }
}
\details{
Dyad independence is determined by checking if all of the constituent parts of the object (formula, ergm terms, etc) are flagged as dyad-independent. 
}
\value{
  \code{TRUE} if the model fit or one implied by the formula is
  dyad-independent; \code{FALSE} otherwise.
}

\keyword{ model }
