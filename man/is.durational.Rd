\name{is.durational}
\alias{is.durational}
\alias{is.durational.formula}
\alias{is.durational.ergm}
\title{
Testing for durational dependent 
}
\description{
These functions test whether an ERGM model or formula is
durational dependent or not.
}
\usage{
\method{is.durational}{character}(object, ...)
\method{is.durational}{ergm}(object, ...)
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
