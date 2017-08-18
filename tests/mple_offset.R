#  File tests/mple_offset.R in package ergm, part of the Statnet suite
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
data(sampson)

set.seed(0)

total.theta <- coef(ergm(samplike~edges))
offset.theta <- pi

stopifnot(all.equal(total.theta, coef(ergm(samplike~edges+offset(edges), offset.coef=c(pi)))[1],
          total.theta-offset.theta, tolerance=0.00001))

stopifnot(all.equal(total.theta, coef(ergm(samplike~offset(edges)+edges, offset.coef=c(pi)))[1], 
          total.theta-offset.theta, tolerance=0.00001))

data(florentine)
boo<-flomarriage
boo[1:3,]<-0
foo <- suppressWarnings(
  ergm(flomarriage~edges+offset(edgecov(boo))+gwesp(0.25,fixed=T),offset.coef=20))
if (max(abs(foo$coef), na.rm=T) > 20) stop("MPLE + offset error")

}, "MPLE + offset")
