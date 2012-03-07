#  File ergm/tests/nodrop.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
### Tests to make sure drop=FALSE works.
library(ergm)
data(sampson)

## Shouldn't need to drop.
# MPLE
summary(ergm(samplike~edges, control=control.ergm(drop=FALSE)))
# MCMC
summary(ergm(samplike~edges, control=control.ergm(drop=FALSE, force.main=TRUE)))

## Empty network.
y0 <- network.initialize(10)
# MPLE
summary(ergm(y0~edges, control=control.ergm(drop=FALSE)))
# MCMC
summary(ergm(y0~edges, control=control.ergm(drop=FALSE, force.main=TRUE, init=0)))

## Full network.
y1 <- as.network(matrix(1,10,10))
# MPLE
summary(ergm(y1~edges, control=control.ergm(drop=FALSE)))
# MCMC
summary(ergm(y1~edges, control=control.ergm(drop=FALSE, force.main=TRUE, init=0)))
