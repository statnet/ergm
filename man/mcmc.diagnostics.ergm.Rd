%  File man/mcmc.diagnostics.ergm.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2017 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{mcmc.diagnostics}
\alias{mcmc.diagnostics.ergm}
\alias{mcmc.diagnostics}
\alias{mcmc.diagnostics.default}
\alias{plot.mcmc.list.ergm}
\title{Conduct MCMC diagnostics on an ergm fit}
\description{
  This function prints diagnistic information and creates simple diagnostic plots for the MCMC sampled statistics produced from a fit.
}
\usage{
\method{mcmc.diagnostics}{ergm} (object, 
              center=TRUE,
              esteq=TRUE,
              vars.per.page=3,
              \dots)
\method{plot}{mcmc.list.ergm}(x, 
                main = NULL, 
                vars.per.page = 3, 
                \dots)
}
\arguments{
  \item{object}{An ergm object.  See documentation for
    \code{\link{ergm}}.
  }
  \item{center}{Logical:  If TRUE, ; center the samples
    on the observed statistics.
  }
  \item{esteq}{Logical:  If TRUE, for statistics corresponding to
    curved ERGM terms, summarize the curved statistics by their
    estimating equation values (evaluated at the MLE of any curved parameters) (i.e.,
    \eqn{\eta'_{I}(\hat{\theta})\cdot g_{I}(y)} for \eqn{I} being indices of
    the canonical parameters in question), 
    rather than the canonical (sufficient) vectors of the curved
    statistics (\eqn{g_{I}(y)}).
  }
  \item{vars.per.page}{Number of rows (one variable per row) 
    per plotting page.  Ignored
    if \code{latticeExtra} package is not installed.
  }
  \item{x}{an mcmc.list object containing the mcmc diagnostic samples}
  \item{main}{
   character, main plot heading title
  }
  \item{\dots}{Additional arguments, to be passed to plotting functions.
  }
}

\details{ A pair of plots are produced for each statistic:a trace of the sampled output statistic values on the left and density estimate for each variable in the MCMC chain on the right.  Diagnostics printed to the console include correlations and convergence diagnostics.
  
  Recent changes in the ergm estimation algorithm mean that these plots can no longer be used to ensure that the mean statistics from the model match the observed network statistics. For that functionality, please use the GOF command: \code{gof(object, GOF=~model)}.
  
  In fact, an ergm output \code{object} contains the matrix of
  statistics from the MCMC run as component \code{$sample}.
  This matrix is actually an object of class \code{mcmc} and
  can be used directly in the \code{coda} package to assess MCMC
  convergence. \emph{Hence all MCMC diagnostic methods available
  in \code{coda} are available directly.} See the examples and
  \url{http://www.mrc-bsu.cam.ac.uk/software/bugs/the-bugs-project-winbugs/coda-readme/}.
  

  More information can be found by looking at the documentation of
  \code{\link{ergm}}.
  
}
\references{
%  Warnes, G.W. (2000).  Multi-Chain and Parallel Algorithms for Markov
%  Chain Monte Carlo. Dissertation, Department of Biostatistics,
%  University of Washington,
%
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
\seealso{\code{\link{ergm}}, \code{network} package,
         \code{coda} package, 
         \code{\link{summary.ergm}}
        }
\examples{
\dontrun{
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
}
\keyword{models}
