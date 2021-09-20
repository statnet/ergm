#  File R/print.network.list.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################
#' @describeIn network.list A [print()] method for network lists.
#' @export
print.network.list <- function(x, stats.print=FALSE, ...) {
  summary.network.list(x, stats.print=stats.print, ...)
  invisible(x)
}

