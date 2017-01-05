#  File tests/mple_target.stats.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
library(statnet.common)
opttest({
library(ergm)
library(coda)
n<-500
base.net <- network.initialize(n=n,directed=FALSE)
norm.stats<-c(.7,.1,.5)
target.stats<-norm.stats*n
print(target.stats)
cat("Structural check:\n")
cat("Mean degree:", norm.stats[1]*2,".\n")
cat("Average degree among nodes with degree 2 or higher:", (2*norm.stats[1]-norm.stats[3])/(1-norm.stats[2]-norm.stats[3]),".\n")

ergm.fit<-ergm(base.net~edges+degree(c(0,1)),target.stats=n*norm.stats)
summary(ergm.fit)
ergm.sim<-simulate(ergm.fit,nsim=1000,statsonly=TRUE)

target.stats.sim<-apply(ergm.sim,2,mean)
print(target.stats.sim)
print(effectiveSize(mcmc(ergm.sim)))
print((target.stats.sim-target.stats)/sqrt(apply(ergm.sim,2,var)/effectiveSize(mcmc(ergm.sim))))
}, "MPLE + target statistics")
