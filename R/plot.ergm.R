#  File R/plot.ergm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
#################################################################################
# The <plot.ergm> function does it plotting via the <mcmc.diagnostics> function.
# This function basically serves as a wrapper
#
# --PARAMETERS--
#   x: an ergm object
#   *: a host of parameters, all of which are ignored; for details see the
#      R documentation for <plot.ergm>
#
# --RETURNED--
#   NULL
# 
###############################################################################

#' Deprecated alias for [mcmc.diagnostics()].
#'
#' Originally, a plotting method for latent positions, before latent
#' space modelling was split out into the `latentnet` package, it now
#' prints the MCMC diagnostics, but it may be repurposed in the
#' future.
#' 
#' @param x an [ergm()] object.
#' @param ... additional arguments, passed to [mcmc.diagnostics()].
#'
#' @seealso [mcmc.diagnostics()]
plot.ergm <- function (x, ...)
{
  .Deprecated("mcmc.diagnostics(x,...)")
    if(is.null(x$sample))
      stop("No plotting method is available for ERGM fits without an MCMC sample.")
    else
      mcmc.diagnostics(x,...)
}
