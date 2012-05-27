#  File ergm/tests/mple_offset.R
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

total.theta <- coef(ergm(samplike~edges))
offset.theta <- pi

stopifnot(all.equal(total.theta, coef(ergm(samplike~edges+offset(edges), offset.coef=c(pi)))[1],
          total.theta-offset.theta, tolerance=0.00001))

stopifnot(all.equal(total.theta, coef(ergm(samplike~offset(edges)+edges, offset.coef=c(pi)))[1], 
          total.theta-offset.theta, tolerance=0.00001))
