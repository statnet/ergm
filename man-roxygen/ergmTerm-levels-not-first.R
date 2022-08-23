#  File man-roxygen/ergmTerm-levels-not-first.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
#' @note To include all attribute values is usually not a good idea, because
#'   the sum of all such statistics equals the number of edges and hence a linear
#'   dependency would arise in any model also including \code{edges}. The default,
#'   \code{levels=-1}, is therefore to omit the first (in lexicographic order)
#'   attribute level. To include all levels, pass either \code{levels=TRUE}
#'   (i.e., keep all levels) or \code{levels=NULL} (i.e., do not filter levels).
