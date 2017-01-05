%  File man/control.gof.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2017 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{control.gof}
\alias{control.gof}
\alias{control.gof.formula}
\alias{control.gof.ergm}
\title{ Auxiliary for Controlling ERGM Goodness-of-Fit Evaluation  }
\description{
  Auxiliary function as user interface for fine-tuning ERGM
  Goodness-of-Fit Evaluation.
}
\usage{
control.gof.formula(nsim=100,
                    MCMC.burnin=10000,
                    MCMC.interval=1000,
                    MCMC.prop.weights="default",
                    MCMC.prop.args=list(),

                    MCMC.init.maxedges=20000,
                    MCMC.packagenames=c(),

                    MCMC.runtime.traceplot=FALSE,          
                    network.output="network",

                    seed=NULL,
                    parallel=0,
                    parallel.type=NULL,
                    parallel.version.check=TRUE)

control.gof.ergm(nsim=100,
                 MCMC.burnin=NULL,
                 MCMC.interval=NULL,
                 MCMC.prop.weights=NULL,
                 MCMC.prop.args=NULL,

                 MCMC.init.maxedges=NULL,
                 MCMC.packagenames=NULL,

                 MCMC.runtime.traceplot=FALSE,
                 network.output="network",

                 seed=NULL,
                 parallel=0,
                 parallel.type=NULL,
                 parallel.version.check=TRUE)
}
\details{
    This function is only used within a call to the \code{\link{gof}} function.
    See the \code{usage} section in \code{\link{gof}} for details.
}
\arguments{
  \item{nsim}{Number of networks to be randomly drawn using Markov
    chain Monte Carlo.  This sample of networks provides the
    basis for comparing the model to the observed network.
  }
  \item{MCMC.burnin}{Number of proposals before any MCMC sampling
    is done. It typically is set to a fairly large number.
  }
  \item{MCMC.interval}{Number of proposals between sampled statistics.
  }
  \item{MCMC.prop.weights}{Specifies the proposal distribution used in the MCMC
    Metropolis-Hastings algorithm.  Possible choices are \code{"TNT"} or 
    \code{"random"}; the \code{"default"} is one of these two, depending on the
    constraints in place (as defined by the \code{constraints}
    argument of the \code{\link{ergm}} function), though not all weights
    may be used with all constraints.  
    The \code{TNT} (tie / no tie) option puts roughly equal weight on selecting a
    dyad with or without a tie as a candidate for toggling, whereas the
    \code{random} option puts equal weight on all possible dyads, though the
    interpretation of \code{random} may change according to the constraints in
    place.  When no constraints are in place, the default is TNT, which 
    appears to improve Markov chain mixing particularly for networks with a
    low edge density, as is typical of many realistic social networks.
  }
  \item{MCMC.prop.args}{An alternative, direct way of specifying additional arguments to proposal.
  }
  \item{MCMC.init.maxedges}{Maximum number of edges expected in network.
  }
  \item{MCMC.packagenames}{Names of packages in which to look for change statistic
    functions in addition to those autodetected. This argument should not be
    needed outside of very strange setups.
  }
  \item{MCMC.runtime.traceplot}{Logical:  If TRUE, plot traceplots of 
    the MCMC sample after every MCMC MLE iteration.
  }
  \item{network.output}{R class with which to output 
    networks. The options are "network" (default) and "edgelist.compressed" 
    (which saves space but only supports networks without vertex attributes)
  }
  \item{seed}{Seed value (integer) for the random number generator.  
    See \code{\link[base]{set.seed}}
  }
  \item{parallel}{Number of threads in which to run the
    sampling. Defaults to 0 (no parallelism). See the entry on
    \link[=ergm-parallel]{parallel processing} for
    details and troubleshooting.
  }
  \item{parallel.type}{API to use for parallel
    processing. Supported values are \code{"MPI"} and
    \code{"PSOCK"}. Defaults to using the \code{parallel} package with PSOCK clusters. See
    \code{\link{ergm-parallel}}
  }
  \item{parallel.version.check}{Logical:  If TRUE, check that the version of
    \code{\link[=ergm-package]{ergm}} running on the slave nodes is the
    same as that running on the master node.
  }
}
\value{
  A list with arguments as components.
}
\seealso{ \code{\link{gof}}. The \code{\link{control.simulate}} 
function performs a 
similar function for
\code{\link{simulate.ergm}}; 
\code{\link{control.ergm}} performs a
similar function for \code{\link{ergm}}.   }
%\keyword{goodness of fit}
\keyword{ models }

