#  File tests/testthat/helper-edge-values.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
add_eattr <- function(nw, attr, f = runif){
  nwn <- substitute(nw)
  nw %e% attr <- f(network.edgecount(nw, na.omit = FALSE))

  eval.parent(call("<-", nwn, nw), 1)
}
