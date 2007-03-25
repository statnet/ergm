# Remove isolated vertices (and associated edges) from the network.
#
delete.isolates<-function(x){
  #Check to be sure we were called with a network
  if(!is.network(x))
    stop("delete.isolates requires an argument of class network.")

  isolates <- (1:network.size(x))[is.isolated(x)]
  if(length(isolates)>0){
    invisible(delete.vertices(x,isolates))
  }else{
    invisible(x)
  }
}
# Remove isolated vertices (and associated edges) from the network.
#
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
central.network<-function(x){
  #Check to be sure we were called with a network
  if(!is.network(x))
    stop("central.network requires an argument of class network.")

# require(sna, quietly=TRUE, warn.conflicts=FALSE)
  xd <- network.copy(x)
  delete.isolates(xd)
# amat <- network(1*(tcrossprod(as.sociomatrix(xd))>0))
  amat <- as.matrix.network(xd,"edgelist")
  isolates <- unique(amat[,2])
  if(length(isolates)>0){delete.vertices(xd,isolates)}
  amat <- as.matrix.network(xd,"edgelist")
  isolates <- unique(amat[,1])
  if(length(isolates)>0){delete.vertices(xd,isolates)}
  delete.isolates(xd)
  invisible(xd)
}
