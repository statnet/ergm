#  File man-roxygen/ergmTerm-form.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################
#' @param form how to aggregate tie values in a valued ERGM: `"sum"`
#'   (the default) for a statistic of the form \eqn{\sum_{i,j} x_{i,j}
#'   y_{i,j}}{sum[i,j] x[i,j]*y[i,j]}, where \eqn{y_{i,j}}{y[i,j]} is
#'   the value of dyad \eqn{(i,j)} and \eqn{x_{i,j}}{x[i,j]} is the
#'   term's covariate associated with it; and `"nonzero"` with the
#'   edge considered to be present if its value is not 0. See
#'   [`ergmTerm`] for more information.
