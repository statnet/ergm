#  File R/InitErgmTerm.bipartite.degree.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################

#NOTE: a number of undocumented terms have been removed from this file
# the terms still exist on the experimental_terms svn branch




###########  InitErgmTerm.b1mindegree  ###################

InitErgmTerm.b1mindegree <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                       varnames = c("d"),
                       vartypes = c("numeric"),
                       defaultvalues = list(NULL),
                       required = c(TRUE))
  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite")
  nb2 <- get.network.attribute(nw, "n") - nb1
  name <- "b1mindegree"
  coef.names <- paste("b1mindeg", a$d, sep="")
  inputs <- a$d
  emptynwstats <- rep(0, length(a$d))
  if (any(a$d==0)) { # alter emptynwstats
    emptynwstats[a$d==0] <- nb1
  }
  list(name=name, coef.names=coef.names, #name and coef.names: required
       inputs = inputs, emptynwstats=emptynwstats,
       dependence = TRUE, minval = 0
       )
}

###########  InitErgmTerm.b2mindegree  ###################

InitErgmTerm.b2mindegree <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                       varnames = c("d"),
                       vartypes = c("numeric"),
                       defaultvalues = list(NULL),
                       required = c(TRUE))
  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite")
  nb2 <- get.network.attribute(nw, "n") - nb1
  name <- "b2mindegree"
  coef.names <- paste("b2mindeg", a$d, sep="")
  inputs <- a$d
  emptynwstats <- rep(0, length(a$d))
  if (any(a$d==0)) { # alter emptynwstats
    emptynwstats[a$d==0] <- nb2
  }
  list(name=name, coef.names=coef.names, #name and coef.names: required
       inputs = inputs, emptynwstats=emptynwstats,
       dependence = TRUE, minval = 0)
}



