%  File ergm/man/control.simulate.stergm.Rd
%  Part of the statnet package, http://statnet.org
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
%  Copyright 2012 the statnet development team
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{control.simulate.stergm}
\alias{control.simulate.stergm}
\alias{control.simulate.network}
\title{ Auxiliary for Controlling Separable Temporal ERGM Simulation }
\description{
  Auxiliary function as user interface for fine-tuning STERGM simulation.
}
\usage{
control.simulate.stergm(MCMC.burnin=1000,
                        MCMC.prop.weights.form="default",MCMC.prop.args.form=NULL,
                        MCMC.prop.weights.diss="default",MCMC.prop.args.diss=NULL,                                  
                        MCMC.init.maxedges=20000,
                        MCMC.packagenames="ergm",
                        MCMC.init.maxchanges=1000000)

control.simulate.stergm(MCMC.burnin=1000,
                        MCMC.prop.weights.form="default",MCMC.prop.args.form=NULL,
                        MCMC.prop.weights.diss="default",MCMC.prop.args.diss=NULL,                                  
                        MCMC.init.maxedges=20000,
                        MCMC.packagenames="ergm",
                        MCMC.init.maxchanges=1000000)

}
\details{
    This function is only used within a call to the \code{\link{simulate}} function.
    See the \code{usage} section in \code{\link{simulate.stergm}} for details.
}
\arguments{
  \item{MCMC.burnin}{Number of Metropolis-Hastings steps per phase
                        (formation and dissolution) per time step.
  }
  \item{MCMC.prop.weights.form, MCMC.prop.weights.diss}{Specifies the proposal distribution used in the MCMC
    Metropolis-Hastings algorithm for formation and dissolution, respectively. Possible choices are \code{"TNT"} or 
    \code{"random"}; the \code{"default"}.  
    The \code{TNT} (tie / no tie) option puts roughly equal weight on selecting a
    dyad with or without a tie as a candidate for toggling, whereas the
    \code{random} option puts equal weight on all possible dyads, though the
    interpretation of \code{random} may change according to the constraints in
    place.  When no constraints are in place, the default is TNT, which 
    appears to improve Markov chain mixing particularly for networks with a
    low edge density, as is typical of many realistic social networks.
  }
  \item{MCMC.prop.args.form, MCMC.prop.args.diss}{An alternative, direct way of specifying additional arguments to proposals. 
  }
  \item{MCMC.init.maxchanges}{Maximum number 
    of toggles changes for which to allocate space.
  }
  \item{MCMC.packagenames}{Names of packages in which change statistic
    functions are found.
  }
  \item{MCMC.init.maxedges}{Maximum number of edges expected in network.
  }
}
\value{
  A list with arguments as components.
}
\seealso{ \code{\link{simulate.stergm}}, \code{\link{simulate.formula}}. 
\code{\link{control.stergm}} performs a 
similar function for
\code{\link{stergm}}.   }
\keyword{models}
