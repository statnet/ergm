#  File tests/u_function.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################

# Private storage: valued
library(ergm)

library(ergm.count)
nw <- network.initialize(4, dir=TRUE)
nw[1,2,names.eval="v",add.edges=TRUE] <- 1
nw[1,3,names.eval="v",add.edges=TRUE] <- 1
nw[3,2,names.eval="v",add.edges=TRUE] <- 1
nw %ergmlhs% "response" <- "v"
s <- summary(nw~sum+test.abs.sum.minus.5+test.abs.sum.minus.5(FALSE)+test.abs.sum.minus.5(FALSE,TRUE)+sociomatrix)
stopifnot(all(abs(s[1]-5)==s[2]) && all(abs(s[1]-5)==s[3]) && all(abs(s[1]-5)==s[4]))
stopifnot(all(s[-(1:4)]==c(as.matrix(nw, attrname="v"))))
sim <- simulate(nw~sum,monitor=~test.abs.sum.minus.5+test.abs.sum.minus.5(FALSE)+test.abs.sum.minus.5(FALSE,TRUE)+sociomatrix, reference=~Poisson, coef=0, nsim=100, control=control.simulate.formula(MCMC.burnin=0,MCMC.interval=1))
s <- attr(sim, "stats")
stopifnot(all(abs(s[,1]-5)==s[,2]) && all(abs(s[,1]-5)==s[,3]) && all(abs(s[,1]-5)==s[,4]))
sim.dyads <- t(sapply(lapply(sim, as.matrix, attrname="v"), c))
stopifnot(all(sim.dyads==s[,-(1:4)]))

