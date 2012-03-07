%  File ergm/man/control.san.Rd
%  Part of the statnet package, http://statnet.org
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
%  Copyright 2012 the statnet development team
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{control.san}
\alias{control.san}
\title{ Auxiliary for Controlling SAN }
\description{
  Auxiliary function as user interface for fine-tuning simulated annealing
  algorithm.
}
\usage{
control.san(coef=NULL,

            SAN.tau=1,
            SAN.invcov=NULL,
            SAN.burnin=10000,
            SAN.interval=10000,
            SAN.init.maxedges=20000,
                      
            SAN.prop.weights="default",
            SAN.prop.args=list(),
            SAN.packagenames="ergm",

            network.output="network",

            seed=NULL,
            parallel=0,
            parallel.type=NULL,
            parallel.version.check=TRUE)
}
\details{
    This function is only used within a call to the \code{\link{san}} function.
    See the \code{usage} section in \code{\link{san}} for details.
}
\arguments{
  \item{coef}{
    Vector of model coefficients used for MCMC simulations, 
    one for each model term.
  }
  \item{SAN.tau}{ Currently unused.
  }
  \item{SAN.invcov}{ 
    Initial inverse covariance matrix used to calculate
    Mahalanobis distance in determining how far a proposed MCMC 
    move is from the \code{target.stats} vector.  
    If NULL, taken to be the covariance
    matrix returned when fitting the MPLE if \code{coef==NULL},
    or the identity matrix otherwise.
  }
  \item{SAN.burnin}{Number of MCMC proposals before any sampling
    is done.
  }
  \item{SAN.interval}{Number of proposals between sampled statistics.
  }
  \item{SAN.init.maxedges}{Maximum number of edges expected in network.
  }
  \item{SAN.prop.weights}{Specifies the method to allocate probabilities of
    being proposed to dyads. Defaults to \code{"default"}, which picks a
    reasonable default for the specified constraint. 
    Other possible values are
    \code{"TNT"}, \code{"random"}, and \code{"nonobserved"}, though not
    all values may be used
    with all possible constraints.
  }
  \item{SAN.prop.args}{An alternative, direct way of specifying additional 
    arguments to proposal.
  }
  \item{SAN.packagenames}{Names of packages in which 
    changestatistics are found. 
    Currently unused as \code{ergm} is presumed.
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
    \code{"SOCK"}. Defaults to using the \code{snow} package
    default.
  }
  \item{parallel.version.check}{Logical:  If TRUE, check that the version of
    \code{\link[=ergm-package]{ergm}} running on the slave nodes is the
    same as that running on the master node.
  }
}
\value{
  A list with arguments as components.
}
\seealso{ \code{\link{san}}
}
\keyword{models}
