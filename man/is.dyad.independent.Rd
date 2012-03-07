%  File ergm/man/is.dyad.independent.Rd
%  Part of the statnet package, http://statnet.org
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
%  Copyright 2012 the statnet development team
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{is.dyad.independent}
\alias{is.dyad.independent}
\alias{is.dyad.independent.formula}
\alias{is.dyad.independent.ergm}
\title{
Testing for dyad-independence
}
\description{
These functions test whether an ERGM fit or formula is dyad-independent.
}
\usage{
\method{is.dyad.independent}{ergm}(object, ...)
\method{is.dyad.independent}{formula}(object, 
                                      basis=NULL,
                                      constraints=~.,
                                      \dots)
}
\arguments{
  \item{object}{
    An \code{\link{ergm}} object or an ERGM formula.
  }
  \item{basis, constraints}{
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
