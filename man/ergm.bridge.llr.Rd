\name{ergm.bridge.llr}
\alias{ergm.bridge.llr}
\title{
A simple implementation of bridge sampling to evaluate
log-likelihood-ratio between two ERGM configurations 
}
\description{
This function uses bridge sampling with geometric spacing to estimate the difference between the log-likelihoods of
two parameter configurations for an ERGM via repeated calls to
\code{\link{simulate.formula.ergm}}.
}
\usage{
ergm.bridge.llr(object, response = NULL, from, to, nsteps = 20,
sample.size = 10000, burnin = 10000, basis = NULL, verbose = FALSE, llronly = FALSE, ...)
}
\arguments{
  \item{object}{
    A model formula. See \code{\link{ergm}} for details.
  }
  \item{response, basis}{
  See \code{\link{simulate.ergm}}.
}
  \item{from}{
  The initial parameter configuration.
}
  \item{to}{
  The final parameter configuration.
}
  \item{nsteps}{
  Number of geometric bridges to use.
}
  \item{sample.size}{
  Total number of MCMC draws to use (to be divided up among the bridges, so
  each bridge gets \code{sample.size/nsteps} draws.
}
  \item{sample.size}{
  Total number of steps draws to discart for each bridge (to be divided up among the bridges, so
  each bridge after the first gets \code{sample.size/sqrt(nsteps)}.
}

\item{llronly}{
  Whether only the estiamted log-ratio should be returned. (Defaults to \code{FALSE}.)
}
\item{\dots}{
  Further arguments to \code{\link{simulate.formula.ergm}}.
  }
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
\code{\link{simulate.formula.ergm}}
}
\keyword{ model }

