#  File R/is.na.ergm.R in package ergm, part of the Statnet suite of packages
#  for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################
#' @describeIn ergm Return `TRUE` if the ERGM was fit to a partially observed network and/or an observational process, such as missing (`NA`) dyads.
#' @export
is.na.ergm <- function(x){
  NVL(x$info$obs, !is.null(x$constrained.obs))
}

#' @describeIn ergm Alias to the `is.na()` method.
#' @export
anyNA.ergm <- function(x, ...){
  is.na(x)
}
