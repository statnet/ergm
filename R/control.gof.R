#  File R/control.gof.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################
#=================================================================
# This file contains the 2 following functions for controlling
# goodness-of-fit computations
#            <control.gof.ergm>
#            <control.gof.formula>
#=================================================================



########################################################################
# Both of the <control.gof.X> functions return a control list for
# customing the fitting procedure used by the gof code
#
# --PARAMETERS--
#   prop.weights  : specifies the method used to allocate probabilities
#                  of being proposed to dyads; options are "TNT",
#                   "random", "nonobserved" and "default"; default=
#                   NULL if X is an ergm (which then uses the weights
#                   that the ergm was fit by); default="default" if
#                   X is a formula (which picks a reasonable default
#                   considering any constraints)
#
# --IGNORED--
#   prop.args     : an alternative, direct way of specifying additional
#                   arguments to proposal; as far as I can tell, the
#                   only use for 'prop.args' is to supply the name
#                   of a nodal attribute for use in the
#                   <InitErgmProposal.nobetweengroupties> function, but this
#                   function will never be called in the path from
#                   <ergm.gof> which is the only code using this
#                   control list.
#   maxchanges    : ??; default=1000000
#
# --RETURNED--
#   a list of the above parameters
#
#########################################################################


#' Auxiliary for Controlling ERGM Goodness-of-Fit Evaluation
#' 
#' Auxiliary function as user interface for fine-tuning ERGM Goodness-of-Fit
#' Evaluation.
#' 
#' This function is only used within a call to the \code{\link{gof}} function.
#' See the \code{usage} section in \code{\link{gof}} for details.
#' 
#' @aliases control.gof control.gof.formula control.gof.ergm
#' @param nsim Number of networks to be randomly drawn using Markov chain Monte
#' Carlo.  This sample of networks provides the basis for comparing the model
#' to the observed network.
#' @param MCMC.burnin Number of proposals before any MCMC sampling is done. It
#' typically is set to a fairly large number.
#' @param MCMC.interval Number of proposals between sampled statistics.
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
#' @param MCMC.init.maxedges Maximum number of edges expected in network.
#' @param MCMC.runtime.traceplot Logical: If TRUE, plot traceplots of the MCMC
#' sample after every MCMC MLE iteration.
#' @param network.output R class with which to output networks. The options are
#' "network" (default) and "edgelist.compressed" (which saves space but only
#' supports networks without vertex attributes)
#' @template control_MCMC_parallel
#' @template seed
#' @template control_MCMC_packagenames
#' @return A list with arguments as components.
#' @seealso \code{\link{gof}}. The \code{\link{control.simulate}} function
#' performs a similar function for \code{\link{simulate.ergm}};
#' \code{\link{control.ergm}} performs a similar function for
#' \code{\link{ergm}}.
#' @name control.gof
#' @export control.gof.ergm
control.gof.formula<-function(nsim=100,
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
                              parallel.version.check=TRUE){
  control<-list()
  for(arg in names(formals(sys.function())))
    control[arg]<-list(get(arg))
  
  set.control.class("control.gof.formula")
}

#' @rdname control.gof
#'
#' @description The `control.gof.ergm` version is intended to be used
#'   with [gof.ergm()] specifically and will "inherit" as many control
#'   parameters from [`ergm`] fit as possible().
#'  
#' @export control.gof.formula
control.gof.ergm<-function(nsim=100,
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
                           parallel.version.check=TRUE){
  control<-list()
  for(arg in names(formals(sys.function())))
    control[arg]<-list(get(arg))

  set.control.class("control.gof.ergm")
}
