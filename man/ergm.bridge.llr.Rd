%  File man/ergm.bridge.llr.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2015 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm.bridge.llr}
\alias{ergm.bridge.llr}
\alias{ergm.bridge.0.llk}
\title{
A simple implementation of bridge sampling to evaluate
log-likelihood-ratio between two ERGM configurations 
}
\description{
\code{ergm.bridge.llr} uses bridge sampling with geometric spacing to 
estimate the difference between the log-likelihoods of
two parameter vectors for an ERGM via repeated calls to
\code{\link{simulate.formula.ergm}}.

\code{ergm.bridge.0.llk} is a convenience wrapper around \code{ergm.bridge.llr}: returns the
log-likelihood of configuration `theta' \emph{relative to the reference
measure}. That is, the configuration with theta=0 is defined as
having log-likelihood of 0

See also \code{\link{ergm.bridge.dindstart.llk}} to use dyad-independent ERGM as a staring point.
}
\usage{
  ergm.bridge.llr(object, 
                  response=NULL, 
                  constraints=~., 
                  from, 
                  to, 
                  basis=NULL, 
                  verbose=FALSE, 
                  \dots, 
                  llronly=FALSE, 
                  control=control.ergm.bridge())
                  
  ergm.bridge.0.llk(object, 
                    response=response, 
                    coef, 
                    ..., 
                    llkonly=TRUE, 
                    control=control.ergm.bridge())
}
\arguments{
  \item{object}{
    A model formula. See \code{\link{ergm}} for details.
  }
  \item{response}{Not for release.
  }
  \item{constraints}{A one-sided formula specifying one or more constraints
    on the support of the distribution of the networks being
    simulated. See the documentation for a similar argument for
    \code{\link{ergm}} for more information. For
    \code{simulate.formula}, defaults to no constraints. For
    \code{simulate.ergm}, defaults to using the same constraints as
    those with which \code{object} was fitted.
  }  
  \item{from, to}{
    The initial and final parameter vectors.
  }
  \item{basis}{An optional \code{\link[network]{network}} object to start
    the Markov chain.  If omitted, the default is the left-hand-side of 
    the \code{object}.
  }
  \item{verbose}{Logical:  If TRUE, print detailed information.
  }
  \item{\dots}{
    Further arguments to \code{ergm.bridge.llr} and \code{\link{simulate.formula.ergm}}.
  }
  \item{llronly}{Logical:  If TRUE, only the estiamted log-ratio will 
    be returned.
  }
  \item{control}{Control arguments.  See 
    \code{\link{control.ergm.bridge}} for details.
  }
  \item{coef}{A vector of coefficients for the configuration of interest.}
  \item{llkonly}{Whether only the estiamted log-likelihood should be returned. (Defaults to TRUE.)}
}
\value{
  If \code{llronly=TRUE}, returns the scalar
  log-likelihood-ratio. Otherwise, returns a list with the following components:
  \item{llr}{The estimated log-ratio.}
  \item{llrs}{The estimated log-ratios for each of the \code{nsteps}
    bridges.}
  \item{path}{A numeric matrix with nsteps rows, with each row being the
    respective bridge's parameter configuration.}
  \item{stats}{A numeric matrix with nsteps rows, with each row being the
    respective bridge's vector of simulated statistics.}
  \item{Dtheta.Du}{The gradient vector of the parameter values with
    respect to position of the bridge.}
}
\references{
Hunter, D. R. and Handcock, M. S. (2006)
    \emph{Inference in curved exponential family models for networks},
   Journal of Computational and Graphical Statistics.

}

\seealso{
\code{\link{simulate.formula.ergm}}, \code{\link{ergm.bridge.dindstart.llk}}
}
\keyword{ model }

