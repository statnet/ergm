#  File tests/constrain_degrees_edges.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################

# Also, while we are at it, test the options setting in .onLoad().
options(ergm.eval.loglik=FALSE) # .onLoad() should not clobber this.

library(ergm)

# Check that options are either set to default or preserved.
stopifnot(getOption("ergm.eval.loglik")==FALSE)
stopifnot(getOption("ergm.loglik.warn_dyads")==TRUE)
