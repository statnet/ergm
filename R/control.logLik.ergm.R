#  File R/control.logLik.ergm.R in package ergm, part of the Statnet suite of
#  packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################


#' @rdname control.ergm.bridge
#'
#' @description By default, the bridge sampler inherits its control
#'   parameters from the [ergm()] fit; `control.logLik.ergm()` allows
#'   the user to selectively override them.
#'
#' @details `control.logLik.ergm()` is only used within a call to the
#'   [logLik.ergm()].
#'
#' @seealso [logLik.ergm()]
#' @export control.logLik.ergm
control.logLik.ergm<-function(bridge.nsteps=16,
                              bridge.target.se=NULL,
                              bridge.bidirectional = TRUE,

                              drop = NULL,

                              MCMC.burnin=NULL,
                              MCMC.interval=NULL,
                              MCMC.samplesize=NULL,
                              obs.MCMC.samplesize=MCMC.samplesize,
                              obs.MCMC.interval=MCMC.interval,
                              obs.MCMC.burnin=MCMC.burnin,
                              
                              MCMC.prop=NULL,
                              MCMC.prop.weights=NULL,
                              MCMC.prop.args=NULL,
                              obs.MCMC.prop=MCMC.prop,
                              obs.MCMC.prop.weights=MCMC.prop.weights,
                              obs.MCMC.prop.args=MCMC.prop.args,

                              MCMC.maxedges=Inf,
                              MCMC.packagenames=NULL,
                              
                              term.options=NULL,
                              seed=NULL,
                              parallel=NULL,
                              parallel.type=NULL,
                              parallel.version.check=TRUE,
                              parallel.inherit.MT=FALSE,

                              ...
){
  old.controls <- list(nsteps = bridge.nsteps)
  control <- handle.controls("control.logLik.ergm", ...)
  set.control.class("control.logLik.ergm")
}
