#  File R/edgelist.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
as.edgelist <- function(x, ...) UseMethod("as.edgelist")

as.edgelist.network <- function(x, attrname = NULL, as.sna.edgelist = FALSE, inverted = NULL, ...){
  as.edgelist(as.matrix.network.edgelist(x, attrname=attrname, as.sna.edgelist=as.sna.edgelist,...), n=network.size(x), directed=is.directed(x), bipartite=if(is.bipartite(x)) x%n%"bipartite" else FALSE, loops=has.loops(x), inverted=NVL(inverted, NVL(x%n%"inverted", FALSE)))
}


as.edgelist.matrix <- function(x, n, directed=TRUE, bipartite=FALSE, loops=FALSE, inverted=FALSE, ...){
  if(!directed) x[,1:2] <- cbind(pmin(x[,1],x[,2]),pmax(x[,1],x[,2]))
  if(!loops) x <- x[x[,1]!=x[,2],,drop=FALSE]
  if(bipartite) x <- x[(x[,1]<=bipartite)!=(x[,2]<=bipartite),,drop=FALSE]
  x <- unique(x[order(x[,1],x[,2]),,drop=FALSE])
  attr(x,"n") <- n
  attr(x,"directed") <- directed
  attr(x,"bipartite") <- bipartite
  attr(x,"loops") <- loops
  attr(x,"inverted") <- inverted
  x
}


