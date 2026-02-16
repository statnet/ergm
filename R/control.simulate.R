#  File R/control.simulate.R in package ergm, part of the Statnet suite of
#  packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################

#' Auxiliary for Controlling ERGM Simulation
#' 
#' Auxiliary function as user interface for fine-tuning ERGM
#' simulation. `control.simulate`, `control.simulate.formula`, and
#' `control.simulate.formula.ergm` are all aliases for the same
#' function.
#' 
#' This function is only used within a call to the ERGM [simulate()]
#' function.  See the Usage section in [simulate.ergm()] for
#' details.
#'
#' @templateVar MCMCType MCMC
#'
#' @template control_MCMC_prop
#' @param MCMC.burnin Number of proposals before any MCMC sampling is done. It
#' typically is set to a fairly large number.
#' @param MCMC.interval Number of proposals between sampled statistics.
#'
#' @template control_MCMC.batch
#'
#' @param MCMC.scale For `control.simulate.ergm()` inheriting
#'   `MCMC.burnin` and `MCMC.interval` from the [`ergm`] fit, the
#'   multiplier for the inherited values. This can be useful because
#'   MCMC parameters used in the fit are tuned to generate a specific
#'   effective sample size for the sufficient statistic in a large
#'   MCMC sample, so the inherited values might not generate
#'   independent realisations.
#'
#' @template control_MCMC_effectiveSize
#' 
#' @template control_MCMC_maxedges
#' @param MCMC.runtime.traceplot Logical: If `TRUE`, plot traceplots of the MCMC
#' sample.
#' @param network.output R class with which to output networks. The options are
#' "network" (default) and "edgelist.compressed" (which saves space but only
#' supports networks without vertex attributes)
#' @template term_options
#' @template control_MCMC_parallel
#' @template control_MCMC_packagenames
#' @template control_dots
#' @return A list with arguments as components.
#' @seealso [simulate.ergm()], [simulate.formula()].
#' [control.ergm()] performs a similar function for
#' [ergm()]; [control.gof()] performs a similar function
#' for [gof()].
#'
#' @name control.simulate.ergm
#' @keywords models
#' @export control.simulate.formula.ergm
control.simulate.formula.ergm<-function(MCMC.burnin=MCMC.interval*16,
                                        MCMC.interval=1024,
                                        MCMC.prop=trim_env(~sparse + .triadic),
                                        MCMC.prop.weights="default",
                                        MCMC.prop.args=list(),

                                        MCMC.batch=NULL,

                                        MCMC.effectiveSize=NULL,
                                        MCMC.effectiveSize.damp=10,
                                        MCMC.effectiveSize.maxruns=1000,
                                        MCMC.effectiveSize.burnin.pval=0.2,
                                        MCMC.effectiveSize.burnin.min=0.05,
                                        MCMC.effectiveSize.burnin.max=0.5,
                                        MCMC.effectiveSize.burnin.nmin=16,
                                        MCMC.effectiveSize.burnin.nmax=128,
                                        MCMC.effectiveSize.burnin.PC=FALSE,
                                        MCMC.effectiveSize.burnin.scl=1024,
                                        MCMC.effectiveSize.order.max=NULL,
                                        
                                        MCMC.maxedges=Inf,
                                        MCMC.packagenames=c(),
                                        
                                        MCMC.runtime.traceplot=FALSE,  
                                        network.output="network",
                                        
                                        term.options=NULL,
                                        
                                        parallel=0,
                                        parallel.type=NULL,
                                        parallel.version.check=TRUE,
                                        parallel.inherit.MT=FALSE,

                                        ...){
  old.controls <- list(
                       maxedges="MCMC.maxedges",
                       prop.weights="MCMC.prop.weights",
                       prop.args="MCMC.prop.args",
                       packagenames="MCMC.packagenames"
                       )

  control <- handle.controls("control.simulate.formula", ...)
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
                                MCMC.scale=1,
                                MCMC.prop=NULL,
                                MCMC.prop.weights=NULL,
                                MCMC.prop.args=NULL,

                                MCMC.batch=NULL,
                                
                                MCMC.effectiveSize=NULL,
                                MCMC.effectiveSize.damp=10,
                                MCMC.effectiveSize.maxruns=1000,
                                MCMC.effectiveSize.burnin.pval=0.2,
                                MCMC.effectiveSize.burnin.min=0.05,
                                MCMC.effectiveSize.burnin.max=0.5,
                                MCMC.effectiveSize.burnin.nmin=16,
                                MCMC.effectiveSize.burnin.nmax=128,
                                MCMC.effectiveSize.burnin.PC=FALSE,
                                MCMC.effectiveSize.burnin.scl=1024,
                                MCMC.effectiveSize.order.max=NULL,
                                
                                MCMC.maxedges=Inf,
                                MCMC.packagenames=NULL,
                                
                                MCMC.runtime.traceplot=FALSE,
                                network.output="network",

                                term.options=NULL,
                                
                                parallel=0,
                                parallel.type=NULL,
                                parallel.version.check=TRUE,
                                parallel.inherit.MT=FALSE,

                                ...){
  old.controls <- list(
                       maxedges="MCMC.maxedges",
                       prop.weights="MCMC.prop.weights",
                       prop.args="MCMC.prop.args",
                       packagenames="MCMC.packagenames"
                       )

  control <- handle.controls("control.simulate.ergm", ...)
  set.control.class("control.simulate.ergm")
}
