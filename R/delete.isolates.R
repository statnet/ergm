#  File R/delete.isolates.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
#================================================================
# This file contains the 3 following functions for converting
# networks into a subgraph of the original graph
#================================================================


##################################################################
# The <delete.isolates> function deletes the isolated nodes from
# a given network.
#
# --PARAMETERS--
#   x: a network
#
# --RETURNED--
#   x: the original network x, with its isolates removed
#
###################################################################

delete.isolates<-function(x){
  #Check to be sure we were called with a network
  if(!is.network(x))
    stop("delete.isolates requires an argument of class network.")

  isolates <- (1:network.size(x))[is.isolate(x)]
  if(length(isolates)>0){
    invisible(delete.vertices(x,isolates))
  }else{
    invisible(x)
  }
}




##################################################################
# The <largest.components> function returns a copy of the given
# network with any components of a size smaller than that specified
# deleted.
#
# --PARAMETERS--
#   x      : a network
#   minsize: the smallest component size that will be kept in x
#
# --RETURNED--
#   xd: the original network x, with the components of size <
#       'minsize' - 1 removed
#
###################################################################

largest.components<-function(x, minsize=4){
  #Check to be sure we were called with a network
  if(!is.network(x))
    stop("largest.components requires an argument of class network.")

  require(sna, quietly=TRUE, warn.conflicts=FALSE)
  xd <- network.copy(x)
  delete.isolates(xd)
  amat <- network(1*(tcrossprod(as.sociomatrix(xd))>0))
  cdist <- component.dist(amat)
# inlarge <- seq(along=cdist$csize)[cdist$csize == max(cdist$csize)]
  inlarge <- seq(along=cdist$csize)[cdist$csize >= minsize]
  isolates <- 1:nrow(amat)
  isolates <- isolates[!(cdist$membership %in% inlarge)]
  if(length(isolates)>0){delete.vertices(xd,isolates)}
  delete.isolates(xd)
  invisible(xd)
}




####################################################################
# The <central.network> function returns an empty graph
#
# --PARAMETERS--
#   x      : a network
#
# --RETURNED--
#   xd: the original network x, with all edges and all nodes removed
#
#####################################################################

central.network<-function(x){
  #Check to be sure we were called with a network
  if(!is.network(x))
    stop("central.network requires an argument of class network.")

# require(sna, quietly=TRUE, warn.conflicts=FALSE)
  xd <- network.copy(x)
  delete.isolates(xd)
# amat <- network(1*(tcrossprod(as.sociomatrix(xd))>0))
  amat <- as.edgelist(xd)
  isolates <- unique(amat[,2])
  if(length(isolates)>0){delete.vertices(xd,isolates)}
  amat <- as.edgelist(xd)
  isolates <- unique(amat[,1])
  if(length(isolates)>0){delete.vertices(xd,isolates)}
  delete.isolates(xd)
  invisible(xd)
}
