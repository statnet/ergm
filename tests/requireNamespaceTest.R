#  File tests/requireNamespaceTest.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################
library(network)
data(flo)

# Also, while we are at it, test the options setting in .onLoad().
options(ergm.eval.loglik=FALSE) # .onLoad() should not clobber this.

requireNamespace('ergm') # Load the namespace, but don't attach the package.

# Check that options are either set to default or preserved.
stopifnot(getOption("ergm.eval.loglik")==FALSE)
stopifnot(getOption("ergm.loglik.warn_dyads")==TRUE)

# run a summary
ergm::summary_formula(as.network(flo)~density)

# try a user-defined ergm term
InitErgmTerm.myedges<-ergm:::InitErgmTerm.edges
ergm::summary_formula(as.network(flo)~myedges)

# actually run ergm()
data(sampson, package="ergm")
fit <- ergm::ergm(samplike~edges)
stopifnot(isTRUE(all.equal(-log(1/(network.edgecount(samplike)/network.dyadcount(samplike))-1), coef(fit), check.attributes=FALSE)))

library(ergm) # Now, attach the package.
# Check that options are either set to default or preserved.
stopifnot(getOption("ergm.eval.loglik")==FALSE)
stopifnot(getOption("ergm.loglik.warn_dyads")==TRUE)
