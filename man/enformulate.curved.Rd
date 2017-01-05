%  File man/enformulate.curved.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2017 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{enformulate.curved}
\alias{enformulate.curved.formula}
\alias{enformulate.curved.ergm}
\alias{enformulate.curved}
\title{
Convert a curved ERGM into a form suitable as initial values for the
same ergm.
}
\description{
The generic \code{enformulate.curved} converts an \code{\link{ergm}} object or
formula of a model with curved terms to the variant in which the curved
parameters embedded into the formula and are removed from the parameter
vector. This is the form required by \code{\link{ergm}} calls.
}
\usage{
\method{enformulate.curved}{ergm}(object, ...)
\method{enformulate.curved}{formula}(object, theta, response=NULL, ...)
}
\arguments{
  \item{object}{
    An \code{\link{ergm}} object or an ERGM formula. The curved terms of
    the given formula (or the formula used in the fit) must have all of their arguments passed by name.
  }
  \item{theta}{
    Curved model parameter configuration.
  }
  \item{response}{Not for release.
  }
  \item{\dots}{
    Unused at this time.
  }
}
\details{
  Because of a current kludge in \code{\link{ergm}}, output from one
  run cannot be directly passed as initial values (\code{control.ergm(init=)}) for
  the next run if any of the terms are curved. One workaround is to
  embed the curved parameters into the formula (while keeping
  \code{fixed=FALSE}) and remove them from \code{control.ergm(init=)}.
  
  This function automates this process for curved ERGM terms included
  with the \code{\link[=ergm-package]{ergm}} package. It does not work with curved terms not
  included in ergm.
}
\value{
  A list with the following components:
  \item{formula}{The formula with curved parameter estimates incorporated.}
  \item{theta}{The coefficient vector with curved parameter estimates removed.}
}

\seealso{
  \code{\link{ergm}},   \code{\link{simulate.ergm}}
}
\examples{
\donttest{
data(sampson)
gest<-ergm(samplike~edges+gwesp(alpha=.5, fixed=FALSE), 
    control=control.ergm(MCMLE.maxit=1))
# Error:
gest2<-try(ergm(gest$formula, control=control.ergm(init=coef(gest), MCMLE.maxit=2)))
print(gest2)

# Works:
tmp<-enformulate.curved(gest)
tmp
gest2<-try(ergm(tmp$formula, control=control.ergm(init=tmp$theta, MCMLE.maxit=2)))
summary(gest2)
}
}
\keyword{ model }
