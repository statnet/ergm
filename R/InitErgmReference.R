#  File R/InitErgmReference.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################
InitErgmReference.Bernoulli <- function(lhs.nw, ...){
  list(name="Bernoulli", init_methods=c("MPLE", "CD", "zeros"))
}


InitErgmReference.StdNormal <- function(lhs.nw, ...){
  list(name="StdNormal", init_methods=c("CD","zeros"))
}

InitErgmReference.Unif <- function(lhs.nw, a, b, ...){
  list(name="Unif", arguments=list(a=a, b=b), init_methods=c("CD","zeros"))
}

InitErgmReference.DiscUnif <- function(lhs.nw, a, b, ...){
  list(name="DiscUnif", arguments=list(a=a, b=b), init_methods=c("CD","zeros"))
}
