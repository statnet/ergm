#  File R/control.simulate.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2018 Statnet Commons
#######################################################################

#' Auxiliary for Controlling ERGM Simulation
#' 
#' Auxiliary function as user interface for fine-tuning ERGM
#' simulation. `control.simulate`, `control.simulate.formula`, and
#' `control.simulate.formula.ergm` are all aliases for the same
#' function.
#' 
#' This function is only used within a call to the \code{\link{simulate}}
#' function.  See the \code{usage} section in \code{\link{simulate.ergm}} for
#' details.
#'
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
#' @param MCMC.burnin Number of proposals before any MCMC sampling is done. It
#' typically is set to a fairly large number.
#' @param MCMC.interval Number of proposals between sampled statistics.
#' @param MCMC.init.maxedges Maximum number of edges expected in network.
#' @param MCMC.runtime.traceplot Logical: If TRUE, plot traceplots of the MCMC
#' sample after every MCMC MLE iteration.
#' @param network.output R class with which to output networks. The options are
#' "network" (default) and "edgelist.compressed" (which saves space but only
#' supports networks without vertex attributes)
#' @template term_options
#' @template control_MCMC_parallel
#' @template control_MCMC_packagenames
#' @param \dots Additional arguments, passed to other functions This argument
#' is helpful because it collects any control parameters that have been
#' deprecated; a warning message is printed in case of deprecated arguments.
#' @return A list with arguments as components.
#' @seealso \code{\link{simulate.ergm}}, \code{\link{simulate.formula}}.
#' \code{\link{control.ergm}} performs a similar function for
#' \code{\link{ergm}}; \code{\link{control.gof}} performs a similar function
#' for \code{\link{gof}}.
#'
#' @name control.simulate.ergm
#' @keywords models
#' @export control.simulate.formula.ergm
control.simulate.formula.ergm<-function(MCMC.burnin=10000,
                                        MCMC.interval=1000,
                                        MCMC.prop.weights="default",
                                        MCMC.prop.args=list(),
                                        
                                        MCMC.init.maxedges=20000,
                                        MCMC.packagenames=c(),
                                        
                                        MCMC.runtime.traceplot=FALSE,  
                                        network.output="network",
                                        
                                        term.options=NULL,
                                        
                                        parallel=0,
                                        parallel.type=NULL,
                                        parallel.version.check=TRUE,
                                        ...){
  old.controls <- list(
                       maxedges="MCMC.init.maxedges",
                       prop.weights="MCMC.prop.weights",
                       prop.args="MCMC.prop.args",
                       packagenames="MCMC.packagenames"
                       )

  control<-list()
  formal.args<-formals(sys.function())
  formal.args[["..."]]<-NULL
  for(arg in names(formal.args))
    control[arg]<-list(get(arg))

  for(arg in names(list(...))){
    if(!is.null(old.controls[[arg]])){
      warning("Passing ",arg," to control.simulate.formula(...) is deprecated and may be removed in a future version. Specify it as control.simulate.formula(",old.controls[[arg]],"=...) instead.")
      control[old.controls[[arg]]]<-list(list(...)[[arg]])
    }else{
      stop("Unrecognized control parameter: ",arg,".")
    }
  }

  set.control.class("control.simulate.formula")
}

#' @rdname control.simulate.ergm
#' @export control.simulate
control.simulate<-control.simulate.formula.ergm
#' @rdname control.simulate.ergm
#' @export control.simulate.formula
control.simulate.formula<-control.simulate.formula.ergm

#' @rdname control.simulate.ergm
#'
#' @description While the others supply a full set of simulation
#'   settings, `control.simulate.ergm` when passed as a control
#'   parameter to [simulate.ergm()] allows some settings to be
#'   inherited from the ERGM stimation while overriding others.
#' 
#' @export control.simulate.ergm
control.simulate.ergm<-function(MCMC.burnin=NULL,
                                MCMC.interval=NULL,
                                MCMC.prop.weights=NULL,
                                MCMC.prop.args=NULL,

                                MCMC.init.maxedges=NULL,
                                MCMC.packagenames=NULL,
                                
                                MCMC.runtime.traceplot=FALSE,
                                network.output="network",

                                term.options=NULL,
                                
                                parallel=0,
                                parallel.type=NULL,
                                parallel.version.check=TRUE,
                                ...){
  old.controls <- list(
                       maxedges="MCMC.init.maxedges",
                       prop.weights="MCMC.prop.weights",
                       prop.args="MCMC.prop.args",
                       packagenames="MCMC.packagenames"
                       )

  control<-list()
  formal.args<-formals(sys.function())
  formal.args[["..."]]<-NULL
  for(arg in names(formal.args))
    control[arg]<-list(get(arg))

  for(arg in names(list(...)))
    if(!is.null(old.controls[[arg]])){
      warning("Passing ",arg," to control.simulate.ergm(...) is deprecated and may be removed in a future version. Specify it as control.simulate.ergm(",old.controls[[arg]],"=...) instead.")
      control[old.controls[[arg]]]<-list(list(...)[[arg]])
    }
 
  set.control.class("control.simulate.ergm")
}
