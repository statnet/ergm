%  File ergm/man/control.logLik.ergm.Rd
%  Part of the statnet package, http://statnet.org
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
%  Copyright 2012 the statnet development team
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{control.logLik.ergm}
\alias{control.logLik.ergm}
\title{ Auxiliary for Controlling logLik.ergm}
\description{
  Auxiliary function as user interface for fine-tuning logLik.ergm
  algorithm, which approximates log likelihood values.
}
\usage{
control.logLik.ergm(nsteps=20,
                    MCMC.burnin=NULL,
                    MCMC.interval=NULL,
                    MCMC.samplesize=NULL,
                    obs.MCMC.samplesize=MCMC.samplesize,
                    obs.MCMC.interval=MCMC.interval,
                    obs.MCMC.burnin=MCMC.burnin,

                    MCMC.prop.weights=NULL,
                    MCMC.prop.args=NULL,

                    MCMC.init.maxedges=NULL,
                    MCMC.packagenames=NULL,

                    seed=NULL)
}
\details{
    This function is only used within a call to the \code{\link{logLik.ergm}} function.
}
\arguments{
  \item{nsteps}{Number of geometric bridges to use.
  }
  \item{MCMC.burnin}{Number of proposals before any MCMC sampling
    is done. It typically is set to a fairly large number.
  }
  \item{MCMC.interval}{Number of proposals between sampled statistics.
  }
  \item{MCMC.samplesize}{Number of network statistics,
    randomly drawn from a given distribution on the set of all networks,
    returned by the Metropolis-Hastings algorithm.
  }
  \item{obs.MCMC.burnin, obs.MCMC.interval, obs.MCMC.samplesize}{
    The \code{obs} versions of these arguments are for the unobserved
    data simulation algorithm.
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
  \item{MCMC.packagenames}{Names of packages in which change statistic
    functions are found.
  }
  \item{seed}{Seed value (integer) for the random number generator.  
    See \code{\link[base]{set.seed}}
  }
}
\value{
  A list with arguments as components.
}
\seealso{ \code{\link{logLik.ergm}}
}
\keyword{models}
