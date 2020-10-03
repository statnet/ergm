#  File tests/fix.curved.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################
library(ergm)
data(sampson)
out<-fix.curved(samplike~edges+gwnsp(decay=.5,fixed=TRUE)+gwesp(decay=.5,fixed=FALSE)+gwodegree(decay=.5,fixed=FALSE)+edges,c(1:7))
stopifnot(out$formula==(samplike ~ edges + gwnsp(decay = 0.5, fixed = TRUE) + gwesp(decay = 4L, fixed = TRUE) + gwodegree(decay = 6L, fixed = TRUE) + edges),
          all(out$theta==c(1,2,3,5,7)))
