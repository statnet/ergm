#  File R/InitErgmReference.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################
InitErgmReference.Bernoulli <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist)
  list(name="Bernoulli", init_methods=c("MPLE", "CD", "zeros"))
}


InitErgmReference.StdNormal <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist)
  list(name="StdNormal", init_methods=c("CD","zeros"))
}

InitErgmReference.Unif <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("a", "b"),
                      vartypes = c("numeric", "numeric"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, TRUE))
  list(name="Unif", arguments=list(a=a$a, b=a$b), init_methods=c("CD","zeros"))
}

InitErgmReference.DiscUnif <- function(nw, arglist, a, b, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("a", "b"),
                      vartypes = c("numeric", "numeric"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, TRUE))
  if(a$a!=round(a$a) || a$b != round(a$b)) ergm_Init_abort(paste("arguments ", sQuote("a"), "and", sQuote("b"), "must be integers"))
  list(name="DiscUnif", arguments=list(a=a$a, b=a$b), init_methods=c("CD","zeros"))
}
