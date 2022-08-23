#  File R/InitErgmTerm.bipartite.degree.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################

#NOTE: a number of undocumented terms have been removed from this file
# the terms still exist on the experimental_terms svn branch




###########  InitErgmTerm.b1mindegree  ###################

#' @templateVar name b1mindegree
#' @title Minimum degree for the first mode in a bipartite network
#' @description This term adds one network statistic to the model for
#'   each element in `d` ; the \eqn{i} th such statistic equals the number of
#'   nodes in the first mode of a bipartite network with at least degree `d[i]` .
#'   The first mode of a bipartite network object is sometimes known as the "actor" mode.
#'   
#' @usage
#' # binary: b1mindegree(d)
#'
#' @param d a vector of distinct integers. 
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-bipartite
#'
#' @concept bipartite
#' @concept undirected
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

#' @templateVar name b2mindegree
#' @title Minimum degree for the second mode in a bipartite network
#' @description This term adds one network statistic to the model for
#'   each element in `d` ; the \eqn{i} th such statistic equals the number of
#'   nodes in the second mode of a bipartite network with at least degree `d[i]` .
#'   The second mode of a bipartite network object is sometimes known as the "event" mode.
#'   
#' @usage
#' # binary: b2mindegree(d)
#'
#' @param d a vector of distinct integers
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-bipartite
#'
#' @concept bipartite
#' @concept undirected
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



