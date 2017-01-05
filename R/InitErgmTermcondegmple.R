#  File R/InitErgmTermcondegmple.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
InitErgmTerm.conddegmple<-function (nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, 
                      varnames = c("type"),
                      vartypes = c("character"),
                      defaultvalues = list("tetrad"),
                      required = c(FALSE))
  ### Process the arguments
  type<-a$type
  ### Construct the list to return
  out=list(name="conddegmple",
       coef.names="conddegmple",
       inputs=NULL,
       emptynwstats=NULL,
       dependence=TRUE)
  out
}
