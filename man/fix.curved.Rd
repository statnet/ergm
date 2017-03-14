%  File man/fix.curved.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2017 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{fix.curved}
\alias{fix.curved.formula}
\alias{fix.curved.ergm}
\alias{fix.curved}
\title{
Convert a curved ERGM into a corresponding ``fixed'' ERGM.
}
\description{
The generic \code{fix.curved} converts an \code{\link{ergm}} object or
formula of a model with curved terms to the variant in which the curved
parameters are fixed. Note that each term has to be treated as a special case.
}
\usage{
\method{fix.curved}{ergm}(object, ...)
\method{fix.curved}{formula}(object, theta, response = NULL, ...)
}
\arguments{
  \item{object}{
    An \code{\link{ergm}} object or an ERGM formula. The curved terms of
    the given formula (or the formula used in the fit) must have all of their arguments passed by name.
}
  \item{theta}{
    Curved model parameter configuration.
}
  \item{response}{
    For valued ERGM, an edge attribute used as the response variable.
}
  \item{\dots}{
    Unused at this time.
}
}
\details{
  Some ERGM terms such as \code{\link{gwesp}} and \code{\link{gwdegree}}
  have two forms: a curved form, for which their decay or similar
  parameters are to be estimated, and whose canonical statistics is a
  vector of the term's components (\code{\link{esp}(1)},
  \code{\link{esp}(2)}, \dots and \code{\link{degree}(1)},
  \code{\link{degree}(2)}, \dots, respectively) and a "fixed" form where
  the decay or similar parameters are fixed, and whose canonical
  statistic is just the term itself. It is often desirable to fit a
  model estimating the curved parameters but simulate the "fixed"
  statistic.

  This function thus takes in a fit or a formula and performs this
  mapping, returning a ``fixed'' model and parameter specification.  It
  only works for curved ERGM terms included with the \code{\link[=ergm-package]{ergm}}
  package. It does not work with curved terms not included in ergm.
}
\value{
  A list with the following components:
  \item{formula}{The ``fixed'' formula.}
  \item{theta}{The ``fixed'' parameter vector.}
}

\seealso{
  \code{\link{ergm}},   \code{\link{simulate.ergm}}
}
\examples{
\donttest{
data(sampson)
gest<-ergm(samplike~edges+gwesp(decay=.5,fixed=FALSE),
           control=control.ergm(MCMLE.maxit=3))
summary(gest)
# A statistic for esp(1),...,esp(16)
simulate(gest,statsonly=TRUE)

tmp<-fix.curved(gest)
tmp
# A gwesp() statistic only
simulate(tmp$formula, coef=tmp$theta, statsonly=TRUE) 
}
}
\keyword{ model }
