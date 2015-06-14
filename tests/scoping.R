#  File tests/scoping.R in package ergm, part of the Statnet suite
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
temp.func <- function(nw,deg) {
  fit <- ergm(nw~edges+degree(deg))
  return(fit)
}
mynet <- network.initialize(50,directed=F)
mynet <- simulate(mynet~edges,coef=-2.5)
print(summary(mynet~edges+degree(0:20)))  
print(temp.func(mynet,3))
mynet <- simulate(mynet~edges,coef=-0.5)
print(summary(mynet~edges+degree(0:20)))  
print(temp.func(mynet,3))
},"scoping issues (Ticket #193)")
