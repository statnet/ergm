#  File R/InitReferences.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
InitReference.Bernoulli <- function(lhs.nw, ...){
  list(name="Bernoulli")  
}

InitReference.Unif <- function(lhs.nw, a, b, ...){
  list(name="Unif", a=a, b=b)  
}

InitReference.DiscUnif <- function(lhs.nw, a, b, ...){
  list(name="DiscUnif", a=a, b=b)
}
