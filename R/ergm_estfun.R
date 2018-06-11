#  File R/ergm_estfun.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2018 Statnet Commons
#######################################################################
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
#' estimating equations (3.1) by Hunter and Handcock (2006) are given
#' instead.
#'
#' @param stats An object representing sample statistics with observed values subtracted out.
#' @param theta Model parameter \eqn{q}-vector.
#' @param model An [`ergm_model`] object or its `etamap` element.
#' @param ... Additional arguments for methods.
#'
#' @return An object of the same class as `stats` containing
#'   \eqn{q}-vectors of estimating function values.
#' @export
ergm.estfun <- function(stats, theta, model, ...){
  UseMethod("ergm.estfun")
}

#' @describeIn ergm.estfun Method for matrices with \eqn{p} columns.
#' @export
ergm.estfun.matrix <- function(stats, theta, model, ...){
  etamap <- if(is(model, "ergm_model")) model$etamap else model
  estf <- t(ergm.etagradmult(theta,t(as.matrix(stats)),etamap))[,!etamap$offsettheta,drop=FALSE]
  colnames(estf) <- (if(is(model, "ergm_model")) param_names(model, FALSE) else names(theta))[!etamap$offsettheta]
  -estf
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
