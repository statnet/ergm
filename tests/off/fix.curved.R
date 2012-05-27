#  File ergm/tests/fix.curved.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
library(ergm)
data(sampson)
out<-fix.curved(samplike~edges+gwnsp(alpha=.5,fixed=TRUE)+gwesp(alpha=.5,fixed=FALSE)+gwodegree(decay=.5,fixed=FALSE)+edges,c(1:7))
stopifnot(out$formula==(samplike ~ edges + gwnsp(alpha = 0.5, fixed = TRUE) + gwesp(alpha = 4L, fixed = TRUE) + gwodegree(decay = 6L, fixed = TRUE) + edges),
          all(out$theta==c(1,2,3,5,7)))
