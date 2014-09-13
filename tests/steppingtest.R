#  File tests/steppingtest.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
library(statnet.common)
opttest({
library(ergm)
library(Rglpk)

####Load the data (provided in the package):
data(ecoli)
form <- ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = T)

m2<-ergm(formula=form, verbose=FALSE, 
        control=control.ergm(main.method="Stepping", Step.MCMC.samplesize=100, Step.gridsize=10000,
        MCMLE.metric="lognormal", MCMC.samplesize=1000, MCMC.burnin=1e+4, MCMC.interval=1000, 
        seed=12345))
if (m2$iterations <5 || m2$iterations > 25) stop("Something fishy in stepping test: Iterations = ", m2$iterations)
}, "Stepping test")
