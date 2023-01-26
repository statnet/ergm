#  File man-roxygen/ergmTerm-x-attrname.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2023 Statnet Commons
################################################################################
#' @param x,attrname either a square matrix of covariates, one for
#'   each possible edge in the network, the name of a network
#'   attribute of covariates, or a network; if the latter, or if the
#'   network attribute named by `x` is itself a network, optional
#'   argument `attrname` provides the name of the quantitative edge
#'   attribute to use for covariate values (in this case, missing
#'   edges in `x` are assigned a covariate value of zero).
