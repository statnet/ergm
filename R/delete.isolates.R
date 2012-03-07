#  File ergm/R/delete.isolates.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
##################################################################
# The <delete.isolates> function deletes the isolated nodes from
# a given network.
##################################################################
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
##################################################################
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
####################################################################
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
