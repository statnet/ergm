#  File ergm/tests/simpletests.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
# Simulate a network with a high number of nodes with outdegree=3 and a low number with indegree=3:
library(ergm)
data(sampson)
m <- simulate.formula(samplike~odegree(3)+idegree(3), coef=c(100,-100))
s <- summary(m~odegree(3)+idegree(3))
if (diff(s) >=0) stop("failed odegree and idegree simulation test.")
 

