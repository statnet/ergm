\name{control.simulate.stergm}
\alias{control.simulate.stergm}
\alias{control.simulate.network}
\title{ Auxiliary for Controlling Separable Temporal ERGM Simulation }
\description{
  Auxiliary function as user interface for fine-tuning STERGM simulation.
}
\usage{
control.simulate.stergm(MCMC.burnin=1000,
                        MCMC.prop.weights="default",
                        MCMC.prop.args=list(),
                        MCMC.init.maxedges=20000,
                        MCMC.packagenames="ergm",

                        summarizestats=FALSE,
                        MCMC.init.maxchanges=1000000)

control.simulate.network(MCMC.burnin=1000,
                         MCMC.prop.weights="default",
                         MCMC.prop.args=list(),
                         MCMC.init.maxedges=20000,
                         MCMC.packagenames="ergm",

                         summarizestats=FALSE,
                         MCMC.init.maxchanges=1000000)
}
\details{
    This function is only used within a call to the \code{\link{simulate}} function.
    See the \code{usage} section in \code{\link{simulate.stergm}} for details.
}
\arguments{
  \item{MCMC.burnin}{Number of proposals before any MCMC sampling
    is done. It typically is set to a fairly large number.
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
  \item{MCMC.init.maxchanges}{Maximum number 
    of changes for which to allocate space.
  }
  \item{MCMC.packagenames}{Names of packages in which change statistic
    functions are found.
  }
  \item{MCMC.init.maxedges}{Maximum number of edges expected in network.
  }
  \item{summarizestats}{logical; Print out a summary of 
    the sufficient statistics of the generated network. This is useful 
    as a diagnostic.
  }
}
\value{
  A list with arguments as components.
}
\seealso{ \code{\link{simulate.stergm}}, \code{\link{simulate.formula}}. 
\code{\link{control.stergm}} performs a 
similar function for
\code{\link{stergm}}.   }
%\keyword{simulation}
\keyword{models}
