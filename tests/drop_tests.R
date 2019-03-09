#  File tests/drop_tests.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################
library(statnet.common)
opttest({
library(ergm)
logit <- function(p) log(p/(1-p))
data(sampson)

# Just one covariate. Note that the .mcmc tests mainly test detection
# and overriding of control$force.main. Note that 1/2 has been subtracted from the "maxed" matrices. This is to test detection of non-0-1 extremeness.

samplike.m <- as.matrix(samplike, matrix.type="adjacency")

maxed.mple <- ergm(samplike~edgecov(samplike.m-1/2))
stopifnot(coef(maxed.mple)==Inf)

maxed.mcmc <- ergm(samplike~edgecov(samplike.m-1/2),control=control.ergm(force.main=TRUE))
stopifnot(coef(maxed.mcmc)==Inf)

mined.mple <- ergm(samplike~edgecov(-samplike.m))
stopifnot(coef(mined.mple)==-Inf)

mined.mcmc <- ergm(samplike~edgecov(-samplike.m),control=control.ergm(force.main=TRUE))
stopifnot(coef(mined.mcmc)==-Inf)

# Now, blank out some of the 1s in the matrix so that you still have a
# dropped term, but now multiple parameters are meaningful.

samplike.m[4:10,4:10] <- 0

truth <- c(logit((network.edgecount(samplike)-sum(samplike.m))/(network.dyadcount(samplike)-sum(samplike.m))),Inf)

maxed.mple <- ergm(samplike~edges+edgecov(samplike.m))
stopifnot(all.equal(truth, coef(maxed.mple),check.attributes=FALSE))

maxed.mcmc <- ergm(samplike~edges+edgecov(samplike.m), control=control.ergm(force.main=TRUE, MCMLE.maxit=10))
stopifnot(all.equal(truth, coef(maxed.mcmc), check.attributes=FALSE,tolerance=0.1))



truth <- c(logit((network.edgecount(samplike)-sum(samplike.m))/(network.dyadcount(samplike)-sum(samplike.m))),-Inf)

mined.mple <- ergm(samplike~edges+edgecov(-samplike.m))
stopifnot(all.equal(truth, coef(mined.mple),check.attributes=FALSE))

mined.mcmc <- ergm(samplike~edges+edgecov(-samplike.m), control=control.ergm(force.main=TRUE, MCMLE.maxit=10))
stopifnot(all.equal(truth, coef(mined.mcmc), check.attributes=FALSE, tolerance=0.1))

# This is mainly to make sure it doesn't crash for dyad-dependent
# and curved terms.
set.seed(1)
y <- network.initialize(10, directed=FALSE)
y[1,2]<-y[2,3]<-y[3,4]<-1
dummy <- ergm(y~edges+triangles+degree(2)+kstar(5)+gwdegree(1,fixed=FALSE),
              control=control.ergm(MCMLE.maxit=3)) # It doesn't seem to stop for a while.
}, "drop")
