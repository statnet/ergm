#  File man-roxygen/ergmTerm-by.R in package ergm, part of the Statnet suite of
#  packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
#' @param by,levels,homophily the optional argument `by` specifies a vertex attribute (see Specifying
#'   Vertex attributes and Levels (`?nodal_attributes`) for details).
#'   If this is specified and `homophily` is `TRUE` ,
#'   then degrees are calculated using the subnetwork consisting of only
#'   edges whose endpoints have the same value of the `by` attribute.
#'   If `by` is specified and
#'   `homophily` is `FALSE` (the default), then separate degree range
#'   statistics are calculated for nodes having each separate
#'   value of the attribute. `levels` selects which levels of by` to include.
