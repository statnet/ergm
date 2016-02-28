#  File tests/samplike_condoutdeg.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
library(ergm)
data(sampson)


degreedist(samplike)

outdegrees <- apply(as.matrix(samplike, m="a"), 1, sum)
table(outdegrees)

efit <- ergm(samplike ~ edges + triangle, estimate="MPLE")
summary(efit)

#
# This fit holds the out degrees fixed
#
efit <- ergm(samplike ~ edges + triangle, constraints=~odegrees,
  control=control.ergm(MCMLE.maxit=3, MCMC.samplesize=10000))
summary(efit)

