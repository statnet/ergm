#  File ergm/R/print.network.list.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
print.network.list <- function(x, stats.print=FALSE, ...) {
  summary.network.list(x, stats.print=stats.print, ...)
  invisible(x)
}

