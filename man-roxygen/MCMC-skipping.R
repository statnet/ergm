#  File man-roxygen/MCMC-skipping.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
#' @section Skipping MCMC iterations (advanced):
#'
#' In some scenarios, it is helpful to forbid certain network
#' configurations from being sampled this can be specified using
#' `constraints`, or by creating an `offset()` term which has value 0
#' if the network is allowed and positive (negative) if the network is
#' not, with offset coefficient set to \code{-\link{Inf}}
#' (\code{+\link{Inf}}). Sometimes, however, a permitted configuration
#' can be reached by "passing through" a forbidden one. For example,
#' if isolates are possible but not nodes with degree exactly 1. Then,
#' an offset term with coefficient [`NaN`] (not [`NA`]!) will cause
#' the MCMC to not terminate as long as the value of that offset term
#' is different from 0.
#'
#' Note that this means that MCMC is not guaranteed to terminate, and
#' there are very few safeguards at this time.
