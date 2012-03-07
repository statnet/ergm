%  File ergm/man/mcmc.diagnostics.ergm.Rd
%  Part of the statnet package, http://statnet.org
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
%  Copyright 2012 the statnet development team
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{mcmc.diagnostics}
\alias{mcmc.diagnostics.ergm}
\alias{mcmc.diagnostics.stergm}
\alias{mcmc.diagnostics}
\title{Conduct MCMC diagnostics on an ergm or stergm fit}
\description{
  This function prints diagnistic information and creates simple diagnostic plots for the MCMC sampled
  statistics produced from a fit.
}
\usage{
\method{mcmc.diagnostics}{ergm} (object, 
                                  center=TRUE,
                                  curved=TRUE,
                                  vars.per.page=3,
                                  \dots) 
\method{mcmc.diagnostics}{stergm} (object,  
                                  center=TRUE,
                                  curved=TRUE,
                                  vars.per.page=3,
                                  \dots) 
}
\arguments{
  \item{object}{An ergm or stergm object.  See documentation for
    \code{\link{ergm}} or \code{\link{stergm}}, respectively.
  }
  \item{center}{Logical:  If TRUE, ; center the samples
    on the observed statistics.
  }
  \item{curved}{Logical:  If TRUE, summarize the curved statistics 
    (evaluated at the MLE of any non-linear parameters), 
    rather than the raw components of the curved statistics.
  }
  \item{vars.per.page}{Number of rows (one variable per row) 
    per plotting page.  Ignored
    if \code{latticeExtra} package is not installed.
  }
  \item{\dots}{Additional arguments, to be passed to plotting functions.
  }
}

\details{ The plots produced are a trace of the sampled output and a
  density estimate for each variable in the chain.  The diagnostics
  printed include correlations and convergence diagnostics.
  
  In fact, an \code{object} contains the matrix of
  statistics from the MCMC run as component \code{$sample}.
  This matrix is actually an object of class \code{mcmc} and
  can be used directly in the \code{coda} package to assess MCMC
  convergence. \emph{Hence all MCMC diagnostic methods available
  in \code{coda} are available directly.} See the examples and
  \url{http://www.mrc-bsu.cam.ac.uk/bugs/classic/coda04/readme.shtml}.

  More information can be found by looking at the documentation of
  \code{\link{ergm}} and \code{\link{stergm}}.
  
}
\references{
  Raftery, A.E. and Lewis, S.M. (1992).  One long run with diagnostics:
  Implementation strategies for Markov chain Monte Carlo. Statistical
  Science, 7, 493-497.

  Raftery, A.E. and Lewis, S.M. (1995).  The number of iterations,
  convergence diagnostics and generic Metropolis algorithms.  In
  Practical Markov Chain Monte Carlo (W.R. Gilks, D.J. Spiegelhalter and
  S. Richardson, eds.). London, U.K.: Chapman and Hall.

  This function is based on the \code{coda} package 
  It is based on the the
  R function \code{raftery.diag} in \code{coda}.  \code{raftery.diag},
  in turn, is based on the FORTRAN program \code{gibbsit} written by
  Steven Lewis which is available from the Statlib archive.
}

\value{\code{\link{mcmc.diagnostics.ergm}} returns 
  some degeneracy information, if it is included in the original
  object.  The function is mainly used for its side effect, which is
  to produce plots and summary output based on those plots.
}
\seealso{\code{\link{ergm}}, \code{\link{stergm}},\code{network} package,
         \code{coda} package, 
         \code{\link{summary.ergm}}
        }
\examples{
#
data(florentine)
#
# test the mcmc.diagnostics function
#
gest <- ergm(flomarriage ~ edges + kstar(2))
summary(gest)

#
# Plot the probabilities first
#
mcmc.diagnostics(gest)
#
# Use coda directly
#
library(coda)
#
plot(gest$sample, ask=FALSE)
#
# A full range of diagnostics is available 
# using codamenu()
#
}
\keyword{models}
