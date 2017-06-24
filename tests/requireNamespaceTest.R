#  File tests/requireNamespaceTest.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
library(network)
data(flo)
requireNamespace('ergm')   #load the namespace, but don't attach the package

# run a summary
ergm::summary.statistics.formula(as.network(flo)~density)

# try a user-defined ergm term
InitErgmTerm.myedges<-ergm:::InitErgmTerm.edges
ergm::summary.statistics.formula(as.network(flo)~myedges)

# try a term that needs to mysteriously access an environment variable
data(sampson, package="ergm")
ergm::summary.statistics.formula(samplike~hammingmix("group"))


# check that we get an appropriate error if no term exists
# should generate an 'unable to locate termed named ... ' error
# ergm::summary.statistics.formula(as.network(flo)~foobar)
