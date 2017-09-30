#  File R/control.logLik.ergm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################


#' Auxiliary for Controlling logLik.ergm
#' 
#' Auxiliary function as user interface for fine-tuning logLik.ergm algorithm,
#' which approximates log likelihood values.
#' 
#' This function is only used within a call to the \code{\link{logLik.ergm}}
#' function.
#' 
#' @param nsteps Number of geometric bridges to use.
#' @param MCMC.burnin Number of proposals before any MCMC sampling is done. It
#' typically is set to a fairly large number.
#' @param MCMC.interval Number of proposals between sampled statistics.
#' @param MCMC.samplesize Number of network statistics, randomly drawn from a
#' given distribution on the set of all networks, returned by the
#' Metropolis-Hastings algorithm.
#' @param obs.MCMC.burnin,obs.MCMC.interval,obs.MCMC.samplesize The \code{obs}
#' versions of these arguments are for the unobserved data simulation
#' algorithm.
#' @param MCMC.prop.weights Specifies the proposal distribution used in the
#' MCMC Metropolis-Hastings algorithm.  Possible choices are \code{"TNT"} or
#' \code{"random"}; the \code{"default"} is one of these two, depending on the
#' constraints in place (as defined by the \code{constraints} argument of the
#' \code{\link{ergm}} function), though not all weights may be used with all
#' constraints.  The \code{TNT} (tie / no tie) option puts roughly equal weight
#' on selecting a dyad with or without a tie as a candidate for toggling,
#' whereas the \code{random} option puts equal weight on all possible dyads,
#' though the interpretation of \code{random} may change according to the
#' constraints in place.  When no constraints are in place, the default is TNT,
#' which appears to improve Markov chain mixing particularly for networks with
#' a low edge density, as is typical of many realistic social networks.
#' @param MCMC.prop.args An alternative, direct way of specifying additional
#' arguments to proposal.
#' @param warn.dyads Whether or not a warning should be issued when sample
#' space constraints render the observed number of dyads ill-defined.
#' @param MCMC.init.maxedges Maximum number of edges expected in network.
#' @param MCMC.packagenames Names of packages in which to look for change
#' statistic functions in addition to those autodetected. This argument should
#' not be needed outside of very strange setups.
#' @param seed Seed value (integer) for the random number generator.  See
#' \code{\link[base]{set.seed}}
#' @param parallel Number of threads in which to run the sampling. Defaults to
#' 0 (no parallelism). See the entry on [parallel processing][ergm-parallel] for details and troubleshooting.
#' @param parallel.type API to use for parallel processing. Supported values
#' are \code{"MPI"} and \code{"PSOCK"}. Defaults to using the \code{parallel}
#' package with PSOCK clusters. See \code{\link{ergm-parallel}}
#' @param parallel.version.check Logical: If TRUE, check that the version of
#' \code{\link[=ergm-package]{ergm}} running on the slave nodes is the same as
#' that running on the master node.
#' @return A list with arguments as components.
#' @seealso \code{\link{logLik.ergm}}
#' @keywords models
#' @export control.logLik.ergm
control.logLik.ergm<-function(nsteps=20,
                              MCMC.burnin=NULL,
                              MCMC.interval=NULL,
                              MCMC.samplesize=NULL,
                              obs.MCMC.samplesize=MCMC.samplesize,
                              obs.MCMC.interval=MCMC.interval,
                              obs.MCMC.burnin=MCMC.burnin,
                              
                              MCMC.prop.weights=NULL,
                              MCMC.prop.args=NULL,

                              warn.dyads=TRUE,

                              MCMC.init.maxedges=NULL,
                              MCMC.packagenames=NULL,
                              
                              seed=NULL,
                              parallel=NULL,
                              parallel.type=NULL,
                              parallel.version.check=TRUE
){

  control<-list()
  formal.args<-formals(sys.function())
  for(arg in names(formal.args))
    control[arg]<-list(get(arg))

  set.control.class("control.logLik.ergm")
}
