%  File man/ergm.bridge.dindstart.llk.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2015 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm.bridge.dindstart.llk}
\alias{ergm.bridge.dindstart.llk}
\title{
  Bridge sampling to estiamte log-likelihood of an ERGM, using a
  dyad-independent ERGM as a staring point.
}
\description{
  This function is a wrapper around \code{\link{ergm.bridge.llr}} that uses a
  dyad-independent ERGM as a starting point for bridge sampling to
  estimate the log-likelihood for a given dyad-dependent model and
  parameter configuration. The dyad-independent model may be specified
  or can be chosen adaptively.
}
\usage{
  ergm.bridge.dindstart.llk(object, 
                            response=NULL,
                            constraints=~.,
                            coef, 
                            dind=NULL, 
                            coef.dind=NULL,  
                            basis=NULL, 
                            \dots, 
                            llkonly=TRUE, 
                            control=control.ergm.bridge())
}
\arguments{
  \item{object}{
    A model formula. See \code{\link{ergm}} for details.
  }
  \item{response}{
    The name of the edge attribute that is the response. Note that it's
    included solely for consistency, since
    \code{ergm.bridge.dindstart.llk} can only handle binary ERGMs.
  }
  \item{constraints}{
    A model constraints formula. See \code{\link{ergm}} for
    details. Note that only constraints that do not induce dyadic
    dependence can be handled by \code{ergm.bridge.dindstart.llk}.
  }
  \item{coef}{
    A vector of coefficients for the configuration of interest.
  }
  \item{dind}{
    A one-sided formula with the dyad-independent model to use as a
    starting point. Defaults to the dyad-independent terms found in the
    formula \code{object} with an overal density term (\code{edges})
    added if not redundant.
  }
  \item{coef.dind}{
    Parameter configuration for the dyad-independent starting
    point. Defaults to the MLE of \code{dind}. 
  }
  \item{basis}{An optional \code{\link[network]{network}} object to start
    the Markov chain.  If omitted, the default is the left-hand-side of 
    the \code{object}.
  }
  \item{\dots}{
    Further arguments to \code{\link{ergm.bridge.llr}} and \code{\link{simulate.formula.ergm}}.
  }
  \item{llkonly}{
    Whether only the estiamted log-likelihood should be returned. (Defaults to \code{TRUE}.)
  }
  \item{control}{Control parameters.  See \code{\link{control.ergm.bridge}}.
  }
}
\value{
  If \code{llkonly=TRUE}, returns the scalar
  log-likelihood. Otherwise, returns a copy of the list returned by
  \code{\link{ergm.bridge.llr}} in addition to the following components:
  \item{llk.dind}{The log-likelihood of the dyad-independence model.}
  \item{llk}{The estimated log-likelihood.}
}
\references{
  Hunter, D. R. and Handcock, M. S. (2006)
  \emph{Inference in curved exponential family models for networks},
  Journal of Computational and Graphical Statistics.
}

\seealso{
  \code{\link{ergm.bridge.llr}}, \code{\link{simulate.formula.ergm}}
}
\keyword{ model }
