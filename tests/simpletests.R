#  File tests/simpletests.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2018 Statnet Commons
#######################################################################
# Simulate a network with a high number of nodes with outdegree=3 and a low number with indegree=3:
library(statnet.common)
opttest({
library(ergm)
data(sampson)
m <- simulate(samplike~odegree(3)+idegree(3), coef=c(100,-100))
s <- summary(m~odegree(3)+idegree(3))
if (diff(s) >=0) stop("failed odegree and idegree simulation test.")
}, "extreme outdegree and indegree simulation test")
