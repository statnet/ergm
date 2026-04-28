#  File R/ergm_estfun.R in package ergm, part of the Statnet suite of packages
#  for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################
#' Compute the Sample Estimating Function Values of an ERGM.
#'
#' The estimating function for an ERGM is the score function: the
#' gradient of the log-likelihood, equalling \eqn{\eta'(\theta)^\top
#' \{g(y)-\mu(\theta)\}}, where \eqn{g(y)} is a \eqn{p}-vector of
#' observed network sufficient statistic, \eqn{\mu(\theta)} is the
#' expected value of the sufficient statistic under the model for
#' parameter value \eqn{\theta}, and \eqn{\eta'(\theta)} is the
#' \eqn{p} by \eqn{q} Jacobian matrix of the mapping from curved
#' parameters to natural parmeters.  If the model is linear, all
#' non-offset statistics are passed. If the model is curved, the score
#' estimating equations (3.1) by \insertCite{HuHa06i;textual}{ergm} are given
#' instead.
#'
#' @param stats An object representing sample statistics with observed values subtracted out.
#' @param theta Model parameter \eqn{q}-vector.
#' @param model An [`ergm_model`] object or its `etamap` element.
#' @param ... Additional arguments for methods.
#'
#' @return An object of the same class as `stats` containing
#'   \eqn{q}-vectors of estimating function values.
#' @keywords internal
#' @export
ergm.estfun <- function(stats, theta, model, ...){
  UseMethod("ergm.estfun")
}

.rightsize_theta_stats <- function(model, stats, theta) {
  etamap <- if (is(model, "ergm_model")) model$etamap else model
  if (any(ot <- etamap$offsettheta)) {
    etamap_no <- deoffset.etamap(etamap, theta)
    if (length(theta) == length(ot)) theta <- theta[!ot]
    if (any(om <- etamap$offsetmap)) {
      stats <-
        if (is.null(ns <- ncol(stats))) {
          if(length(stats) == length(om)) stats[!om] else stats
        } else {
          as.matrix(if (ns == length(om)) stats[, !om, drop = FALSE] else stats)
        }
    }
  } else etamap_no <- etamap

  if (is(model, "ergm_model")) names(theta) <- param_names(model, FALSE, FALSE)

  assign("model", etamap_no, parent.frame())
  assign("stats", stats, parent.frame())
  assign("theta", theta, parent.frame())

  NULL
}

#' @describeIn ergm.estfun Method for numeric vectors of length \eqn{p}.
#' @export
ergm.estfun.numeric <- function(stats, theta, model, ...) {
  .rightsize_theta_stats(model, stats, theta)

  setNames(-c(ergm.etagradmult(theta, stats, model)), names(theta))
}

#' @describeIn ergm.estfun Method for matrices with \eqn{p} columns.
#' @export
ergm.estfun.matrix <- function(stats, theta, model, ...){
  .rightsize_theta_stats(model, stats, theta)

  structure(-ergm.etagradmultt(theta, stats, model),
            dimnames = list(rownames(stats), names(theta)))
}

#' @describeIn ergm.estfun Method for [`mcmc`] objects with \eqn{p} variables.
#' @export
ergm.estfun.mcmc <- function(stats, theta, model, ...){
  mcmc(ergm.estfun(as.matrix(stats), theta, model), start=start(stats), end=end(stats), thin=thin(stats))
}

#' @describeIn ergm.estfun Method for  [`mcmc.list`] objects with \eqn{p} variables.
#' @export
ergm.estfun.mcmc.list <- function(stats, theta, model, ...){
  lapply.mcmc.list(stats, ergm.estfun, theta, model)
}
