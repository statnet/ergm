\name{ergm.bridge.dindstart.llk}
\alias{ergm.bridge.dindstart.llk}
\title{
  Bridge sampling to estiamte log-likelihood of an ERGM, using a
  dyad-independent ERGM as a staring point.
}
\description{
  This function is a wrapper around \code{\link{bergm.bridge.llr}} that uses a
  dyad-independent ERGM as a starting point for bridge sampling to
  estimate the log-likelihood for a given dyad-dependent model and
  parameter configuration. The dyad-independent model may be specified
  or can be chosen adaptively.
}
\usage{
ergm.bridge.dindstart.llk(object, response = NULL, theta, nsteps, dind = NULL, theta.dind = NULL, basis = NULL, llkonly = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    A model formula. See \code{\link{ergm}} for details.
    }
  \item{response, basis, \dots}{
    
  See \code{\link{simulate.ergm}}.
}
  \item{theta}{
    A vector of coefficients for the configuration of interest.
}
  \item{nsteps}{
    Number of geometric bridges to use.
}
  \item{dind}{
    A one-sided formula with the dyad-independent model to use as a
    starting point. Defaults to the dyad-independent terms found in the
    formula \code{object} with an overal density term (\code{edges})
    added if not redundant.
}
  \item{theta.dind}{
    Parameter configuration for the dyad-independent starting
    point. Defaults to the MLE of \code{dind}. 
}
  \item{llkonly}{
      Whether only the estiamted log-likelihood should be returned. (Defaults to \code{TRUE}.)
    }
    \item{\dots}{
      Further arguments to \code{\link{ergm.bridge.llr}} and \code{\link{simulate.formula.ergm}}.
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
\author{
Pavel N. Krivitsky
}
\seealso{
\code{\link{ergm.bridge.llr}}, \code{\link{simulate.formula.ergm}}
}
\keyword{ model }
