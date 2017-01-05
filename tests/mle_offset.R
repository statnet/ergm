#  File tests/mle_offset.R in package ergm, part of the Statnet suite
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
data(florentine)
 
fit1 <- ergm(flomarriage ~ offset(edges) + kstar(2:3), offset.coef=-1, control=control.ergm(MCMLE.maxit=3))

print(summary(fit1))

fit2 <- ergm(flomarriage ~ edges + offset(kstar(2:3)), offset.coef=c(1,-1), control=control.ergm(MCMLE.maxit=3))

print(summary(fit2))
}, "MLE + offset")
