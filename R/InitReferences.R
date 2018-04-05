#  File R/InitErgmReferences.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
InitErgmReference.Bernoulli <- function(lhs.nw, ...){
  list(name="Bernoulli")  
}


InitErgmReference.StdNormal <- function(lhs.nw, ...){
  list(name="StdNormal")  
}

InitErgmReference.Unif <- function(lhs.nw, a, b, ...){
  list(name="Unif", a=a, b=b)  
}

InitErgmReference.DiscUnif <- function(lhs.nw, a, b, ...){
  list(name="DiscUnif", a=a, b=b)
}
