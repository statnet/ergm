#  File R/control.san.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################


#' Auxiliary for Controlling SAN
#' 
#' Auxiliary function as user interface for fine-tuning simulated annealing
#' algorithm.
#' 
#' This function is only used within a call to the \code{\link{san}} function.
#' See the \code{usage} section in \code{\link{san}} for details.
#'
#' @templateVar MCMCType SAN
#'
#' @param SAN.maxit Number of temperature levels to use.
#' 
#' @param SAN.tau Tuning parameter, specifying the temperature of the
#'   process during the *penultimate* iteration. (During the last
#'   iteration, the temperature is set to 0, resulting in a greedy
#'   search, and during the previous iterations, the temperature is
#'   set to `SAN.tau*(iterations left after this one)`.
#' 
#' @param SAN.invcov Initial inverse covariance matrix used to
#'   calculate Mahalanobis distance in determining how far a proposed
#'   MCMC move is from the \code{target.stats} vector.  If `NULL`,
#'   initially set to the identity matrix. In either case, during
#'   subsequent runs, it is estimated empirically.
#'
#' @param SAN.invcov.diag Whether to only use the diagonal of the
#'   covariance matrix. It seems to work better in practice.
#'
#' @param SAN.nsteps.alloc Either a numeric vector or a function of
#'   the number of runs giving a sequence of relative lengths of
#'   simulated annealing runs.
#'
#' @param SAN.nsteps Number of MCMC proposals for all the annealing runs combined.
#' @param SAN.samplesize Number of realisations' statistics to obtain for tuning purposes.
#' @template control_MCMC_prop
#' @param SAN.packagenames Names of packages in which to look for change
#' statistic functions in addition to those autodetected. This argument should
#' not be needed outside of very strange setups.
#' @param SAN.ignore.finite.offsets Whether SAN should ignore (treat as 0) finite offsets.
#' @template term_options
#' @template control_MCMC_parallel
#' @template seed
#' @return A list with arguments as components.
#' @seealso \code{\link{san}}
#' @keywords models
#' @export control.san
control.san<-function(SAN.maxit=4,
                      SAN.tau=1,
                      SAN.invcov=NULL,
                      SAN.invcov.diag=FALSE,
                      SAN.nsteps.alloc=function(nsim) 2^seq_len(nsim),
                      SAN.nsteps=2^19,
                      SAN.samplesize=2^12,
                      
                      SAN.prop=trim_env(~sparse),
                      SAN.prop.weights="default",
                      SAN.prop.args=list(),
                      SAN.packagenames=c(),
                      
                      SAN.ignore.finite.offsets=TRUE,
                      
                      term.options=list(),

                      seed=NULL,
                      parallel=0,
                      parallel.type=NULL,
                      parallel.version.check=TRUE,
                      parallel.inherit.MT=FALSE){

  control <- handle.controls("control.san")
  set.control.class("control.san")
}
