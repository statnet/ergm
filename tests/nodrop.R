#  File tests/nodrop.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
### Tests to make sure drop=FALSE works.
library(statnet.common)
opttest({
library(ergm)
data(sampson)

## Shouldn't need to drop.
# MPLE
summary(ergm(samplike~edges, control=control.ergm(drop=FALSE)))
# MCMC
summary(ergm(samplike~edges, control=control.ergm(drop=FALSE, force.main=TRUE, MCMLE.maxit=10)))

## Empty network.
y0 <- network.initialize(10)
# MPLE
summary(ergm(y0~edges, control=control.ergm(drop=FALSE)))
# MCMC
summary(ergm(y0~edges, control=control.ergm(drop=FALSE, force.main=TRUE, init=0, MCMLE.maxit=10)))

## Full network.
y1 <- as.network(matrix(1,10,10))
# MPLE
summary(ergm(y1~edges, control=control.ergm(drop=FALSE)))
# MCMC
summary(ergm(y1~edges, control=control.ergm(drop=FALSE, force.main=TRUE, init=0, MCMLE.maxit=10)))
}, "drop disabled")
