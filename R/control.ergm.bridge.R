#  File R/control.ergm.bridge.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################


#' Auxiliaries for Controlling [ergm.bridge.llr()] and [logLik.ergm()]
#' 
#' Auxiliary functions as user interfaces for fine-tuning the
#' [ergm.bridge.llr()] algorithm, which approximates log likelihood
#' ratios using bridge sampling.
#' 
#' `control.ergm.bridge()` is only used within a call to the
#' [ergm.bridge.llr()], [ergm.bridge.dindstart.llk()], or
#' [ergm.bridge.0.llk()] functions.
#'
#' @templateVar MCMCType MCMC
#'
#' @param bridge.nsteps Number of geometric bridges to use.
#' @param bridge.target.se If not `NULL`, if the estimated MCMC standard error of the likelihood estimate exceeds this, repeat the bridge sampling, accumulating samples.
#' @param bridge.bidirectional Whether the bridge sampler first bridges from `from` to `to`, then from `to` to `from` (skipping the first burn-in), etc. if multiple attempts are required.
#' @param MCMC.burnin Number of proposals before any MCMC sampling is done. It
#' typically is set to a fairly large number.
#' @param MCMC.burnin.between Number of proposals between the bridges; typically, less and less is needed as the number of steps decreases.
#' @param MCMC.interval Number of proposals between sampled statistics.
#' @param MCMC.samplesize Number of network statistics, randomly drawn from a
#' given distribution on the set of all networks, returned by the
#' Metropolis-Hastings algorithm.
#' @param obs.MCMC.burnin,obs.MCMC.burnin.between,obs.MCMC.interval,obs.MCMC.samplesize The \code{obs}
#' versions of these arguments are for the unobserved data simulation
#' algorithm.
#' @template control_MCMC_prop
#' @param obs.MCMC.prop,obs.MCMC.prop.weights,obs.MCMC.prop.args The `obs` versions of these arguments are for the unobserved data simulation algorithm.
#' @template control_MCMC_maxedges
#' @template term_options
#' @template control_MCMC_parallel
#' @template seed
#' @template control_MCMC_packagenames
#' @template control_dots
#' @return A list with arguments as components.
#' @seealso [ergm.bridge.llr()]
#' @keywords models
#' @export control.ergm.bridge
control.ergm.bridge<-function(bridge.nsteps=16, # Number of geometric bridges to use
                              bridge.target.se=NULL,
                              bridge.bidirectional = TRUE,

                              MCMC.burnin=MCMC.interval*128,
                              MCMC.burnin.between=max(ceiling(MCMC.burnin/sqrt(bridge.nsteps)), MCMC.interval*16),
                              MCMC.interval=128,
                              MCMC.samplesize=16384, # Total number of MCMC draws to use (to be divided up among the bridges, so each bridge gets \code{sample.size/bridge.nsteps} draws.

                              obs.MCMC.burnin=obs.MCMC.interval*128,
                              obs.MCMC.burnin.between=max(ceiling(obs.MCMC.burnin/sqrt(bridge.nsteps)), obs.MCMC.interval*16),
                              obs.MCMC.interval=MCMC.interval,
                              obs.MCMC.samplesize=MCMC.samplesize,

                              MCMC.prop=trim_env(~sparse),
                              MCMC.prop.weights="default",
                              MCMC.prop.args=list(),
                              obs.MCMC.prop=MCMC.prop,
                              obs.MCMC.prop.weights=MCMC.prop.weights,
                              obs.MCMC.prop.args=MCMC.prop.args,

                              MCMC.maxedges=Inf,
                              MCMC.packagenames=c(),
                              
                              term.options=list(),
                              seed=NULL,
                              parallel=0,
                              parallel.type=NULL,
                              parallel.version.check=TRUE,
                              parallel.inherit.MT=FALSE,

                              ...
){
  old.controls <- list(nsteps = bridge.nsteps)

  control <- handle.controls("control.ergm.bridge", ...)
  set.control.class("control.ergm.bridge")
}
