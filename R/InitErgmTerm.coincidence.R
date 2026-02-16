#  File R/InitErgmTerm.coincidence.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################

#' @templateVar name coincidence
#' @title Coincident node count for the second mode in a bipartite (aka two-mode) network
#' @description By default this term adds one
#'   network statistic to the model for each pair of nodes of mode two. It is
#'   equal to the number of (first mode) mutual partners of that pair.
#'   The first mode of a bipartite
#'   network object is sometimes known as the "actor" mode and the seconds as the "event" mode. So this is the number of actors going to both events in the pair. This term can only be
#'   used with undirected bipartite networks.
#'   
#' @usage
#' # binary: coincidence(levels=NULL,active=0)
#'
#' @templateVar explain specifies which pairs of nodes in mode two to include.
#' @template ergmTerm-levels-doco
#' @param active selects pairs for which the observed count is at least `active` . Ignored if `levels` is
#'   specified. (Thus, indices passed as `levels` should correspond to indices when `levels` = NULL and `active` = 0.) 
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-args-3.9.4
#'
#' @concept bipartite
#' @concept undirected
InitErgmTerm.coincidence<-function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                        varnames = c("d","active"),
                        vartypes = c("matrix","numeric"),
                        defaultvalues = list(NULL,0),
                        required = c(FALSE,FALSE))
	levels <- if(!is.null(a$d)) I(transpose(data.frame(node1 = a$d[,1], node2 = a$d[,2]))) else NULL
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                        varnames = c("levels","active"),
                        vartypes = c(ERGM_LEVELS_SPEC,"numeric"),
                        defaultvalues = list(NULL,0),
                        required = c(FALSE,FALSE))
	levels <- a$levels
  }
  nb1 <- b1.size(nw)
  nb2 <- b2.size(nw)

  # coincidence was originally implemented in such a way that the second column varies fastest
  # this behavior is preserved here, although it violates the convention used by other terms
  levels2.grid <- expand.grid(node2 = seq_len(nb2), node1 = seq_len(nb2))
  levels2.grid <- levels2.grid[,c(2,1)] # make the second column vary fastest
  levels2.grid <- levels2.grid[levels2.grid$node1 < levels2.grid$node2,]

  if(!is.null(levels)) {
	levels2.list <- transpose(levels2.grid)
    levels2.sel <- ergm_attr_levels(levels, list(node1 = seq_len(nb2), node2 = seq_len(nb2)), nw, levels2.list)
	
	active <- !is.na(match(levels2.list, levels2.sel))
  } else if(a$active > 0){
   active=summary(nw ~ coincidence(active=0)) > a$active
  }
  if(!is.null(levels) || (a$active > 0)){
   coef.names <- paste("coincidence.",levels2.grid[active,1],".", levels2.grid[active,2], sep="")
   b <- cumsum(active)
   b[active==0] <- 0
  }else{
   b <- seq_len(nrow(levels2.grid))
   coef.names <- paste("coincidence.",levels2.grid[,1],".", levels2.grid[,2], sep="")
  }
  name <- "coincidence"
  inputs <- c(b)
  list(name=name, coef.names=coef.names, inputs=inputs, dependence=TRUE)
}
