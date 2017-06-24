#  File tests/termTests.b12factor.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
library(ergm)



bipnet<-network.initialize(4,bipartite=2,directed=FALSE)
add.edges(bipnet,tail=c(1),head=c(4))
# set an attribute on the the vertices of the second mode
set.vertex.attribute(bipnet,'felines',c('cat','tiger'),v=3:4)

tryCatch(summary(bipnet~b1factor('felines')),error=function(e)cat("expected error message \n"))


if(summary(bipnet~b2factor('felines')) !=1)
	stop("b2factor error A")


