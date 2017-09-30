#  File R/control.san.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
#' Auxiliary for Controlling SAN
#' 
#' Auxiliary function as user interface for fine-tuning simulated annealing
#' algorithm.
#' 
#' This function is only used within a call to the \code{\link{san}} function.
#' See the \code{usage} section in \code{\link{san}} for details.
#' 
#' @param coef Vector of model coefficients used for MCMC simulations, one for
#' each model term.
#' @param SAN.tau Currently unused.
#' @param SAN.invcov Initial inverse covariance matrix used to calculate
#' Mahalanobis distance in determining how far a proposed MCMC move is from the
#' \code{target.stats} vector.  If NULL, taken to be the covariance matrix
#' returned when fitting the MPLE if \code{coef==NULL}, or the identity matrix
#' otherwise.
#' @param SAN.burnin Number of MCMC proposals before any sampling is done.
#' @param SAN.interval Number of proposals between sampled statistics.
#' @param SAN.init.maxedges Maximum number of edges expected in network.
#' @param SAN.prop.weights Specifies the method to allocate probabilities of
#' being proposed to dyads. Defaults to \code{"default"}, which picks a
#' reasonable default for the specified constraint.  Other possible values are
#' \code{"TNT"}, \code{"random"}, and \code{"nonobserved"}, though not all
#' values may be used with all possible constraints.
#' @param SAN.prop.args An alternative, direct way of specifying additional
#' arguments to proposal.
#' @param SAN.packagenames Names of packages in which to look for change
#' statistic functions in addition to those autodetected. This argument should
#' not be needed outside of very strange setups.
#' @param MPLE.max.dyad.types Maximum number of unique values of change
#' statistic vectors, which are the predictors in a logistic regression used to
#' calculate the MPLE.  This calculation uses a compression algorithm that
#' allocates space based on \code{MPLE.max.dyad.types}
#' @param MPLE.samplesize Not currently documented; used in
#' conditional-on-degree version of MPLE.
#' @param network.output R class with which to output networks. The options are
#' "network" (default) and "edgelist.compressed" (which saves space but only
#' supports networks without vertex attributes)
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
#' @seealso \code{\link{san}}
#' @keywords models
#' @export control.san
control.san<-function(coef=NULL,

                      SAN.tau=1,
                      SAN.invcov=NULL,
                      SAN.burnin=100000,
                      SAN.interval=10000,
                      SAN.init.maxedges=20000,
                      
                      SAN.prop.weights="default",
                      SAN.prop.args=list(),
                      SAN.packagenames=c(),
                      
                      MPLE.max.dyad.types=1e6,
                      MPLE.samplesize=50000,

                      network.output="network",

                      seed=NULL,
                      parallel=0,
                      parallel.type=NULL,
                      parallel.version.check=TRUE){
  control<-list()
  for(arg in names(formals(sys.function())))
    control[arg]<-list(get(arg))

  set.control.class("control.san")
}
