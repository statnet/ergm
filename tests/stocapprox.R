#  File tests/simpletests.R in package ergm, part of the Statnet suite
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
  set.seed(2)
  
  data(florentine)
  
  mod.sa = ergm(flomarriage~edges+triangle,control=control.ergm(main.method="Stochastic-Approximation", SA.nsubphases = 6))
  summary(mod.sa)
  
  mod.mcmle = ergm(flomarriage~edges+triangle)
  summary(mod.sa)
  
  stopifnot(all(abs(mod.sa$coef - mod.mcmle$coef) < 0.5))
  
}, "extreme outdegree and indegree simulation test")
