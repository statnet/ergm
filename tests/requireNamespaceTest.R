#  File tests/requireNamespaceTest.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################
library(network)
data(flo)
requireNamespace('ergm')   #load the namespace, but don't attach the package

# run a summary
ergm::summary_formula(as.network(flo)~density)

# try a user-defined ergm term
InitErgmTerm.myedges<-ergm:::InitErgmTerm.edges
ergm::summary_formula(as.network(flo)~myedges)

# actually run ergm()
data(sampson, package="ergm")
fit <- ergm::ergm(samplike~edges)
stopifnot(isTRUE(all.equal(-log(1/(network.edgecount(samplike)/network.dyadcount(samplike))-1), coef(fit), check.attributes=FALSE)))

# check that we get an appropriate error if no term exists
# should generate an 'unable to locate termed named ... ' error
# ergm::summary_formula.formula(as.network(flo)~foobar)
