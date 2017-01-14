#  File tests/fix.curved.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
library(ergm)
data(sampson)
gest<-ergm(samplike~edges+gwidegree(), 
    control=control.ergm(MCMLE.maxit=1, init=c(-0.91624080, 2.19474452, 0.09277566)),eval.loglik=FALSE, verbose=TRUE)
