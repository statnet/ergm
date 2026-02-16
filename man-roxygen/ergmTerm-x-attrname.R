#  File man-roxygen/ergmTerm-x-attrname.R in package ergm, part of the Statnet
#  suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################
#' @param x,attrname a specification for the dyadic covariate: either
#'   one of the following, or the name of a network attribute
#'   containing one of the following: \describe{
#'
#'   \item{a covariate matrix}{with dimensions \eqn{n \times n}{n*n}
#'     for unipartite networks and \eqn{b \times (n-b)}{b*(n-b)} for
#'     bipartite networks; `attrname`, if given, is used to construct
#'     the term name.}
#'
#'   \item{a network object}{with the same size and bipartitedness as
#'     LHS; `attrname`, if given, provides the name of the
#'     quantitative edge attribute to use for covariate values (in
#'     this case, missing edges in `x` are assigned a covariate value
#'     of zero).}
#'
#' }
