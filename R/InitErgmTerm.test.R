#  File R/InitErgmTerm.test.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################
InitErgmTerm..edges_times<-function(nw, arglist, ..., times) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  
  list(name="_edges_times", coef.names=paste("edges_times",times), dependence=FALSE,
       minval = min(times*network.dyadcount(nw,FALSE),0), maxval = max(times*network.dyadcount(nw,FALSE),0), conflicts.constraints="edges", inputs=c(times))
}
