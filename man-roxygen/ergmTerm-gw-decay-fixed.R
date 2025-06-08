#  File man-roxygen/ergmTerm-gw-decay-fixed.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
#' @param decay nonnegative decay parameter for the <%= multiplicand %>; required if `fixed=TRUE` and ignored with a warning otherwise.
#' @param fixed optional argument indicating
#'   whether the `decay` parameter is fixed at the given value, or is to be fit as a curved
#'   exponential-family model (see Hunter and Handcock, 2006). The
#'   default is `FALSE` , which means the scale parameter is not
#'   fixed and thus the model is a curved exponential family. 
