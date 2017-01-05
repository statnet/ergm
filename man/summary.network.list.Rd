%  File man/summary.network.list.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2017 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{summary.network.list}
\alias{summary.network.list}
\alias{print.network.list}
\alias{network.list} 
\title{Summarizing network.list objects}
\usage{
  \method{summary}{network.list}(object, 
                               stats.print=TRUE, 
                               net.print=FALSE, 
                               net.summary=FALSE,
                               \dots)
               
  \method{print}{network.list}(x, stats.print=FALSE, \dots)
}
\arguments{
  \item{object, x}{an object of class \code{network.list}, such 
    as the output from \code{\link{simulate.ergm}}
  }
  \item{stats.print}{Logical:  If TRUE,
    print network statistics.    
  }
  \item{net.print}{Logical:  If TRUE,
    print network overviews.
  }
  \item{net.summary}{Logical:  If TRUE,
    print network summaries.
  }
  \item{\dots}{Additional arguments to be passed to lower-level
    functions.    
  }
}
\description{
\code{\link[base]{summary}} and
\code{\link[base]{print}} methods for class \code{network.list}.
}
\value{
  The \code{summary.network.list} function returns a 
  \code{\link[network]{summary.network}} object.
  The \code{print.summary.list} function calls the
  \code{summary.network.list} function but returns the
  \code{network.list} object.
}
\seealso{ \code{\link{simulate.ergm}}
}
\examples{
# Draw from a Bernoulli model with 16 nodes
# and tie probability 0.1
#
g.use <- network(16, density=0.1, directed=FALSE)
#
# Starting from this network let's draw 3 realizations
# of a model with edges and 2-star terms
#
g.sim <- simulate(~edges+kstar(2), nsim=3, coef=c(-1.8, 0.03),
               basis=g.use, control=control.simulate(
                 MCMC.burnin=100000,
                 MCMC.interval=1000))
print(g.sim)
summary(g.sim)
}
\keyword{regression}
\keyword{models}
