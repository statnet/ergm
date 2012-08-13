as.edgelist <- function(x, ...) UseMethod("as.edgelist")

as.edgelist.network <- function(x, attrname = NULL, as.sna.edgelist = FALSE,...){
  as.edgelist(as.matrix.network.edgelist(x, attrname=attrname, as.sna.edgelist=as.sna.edgelist,...))
}

as.edgelist.matrix <- function(x, directed=TRUE, ...){
  if(!directed) x[,1:2] <- cbind(pmin(x[,1],x[,2]),pmax(x[,1],x[,2]))
  x[order(x[,1],x[,2]),,drop=FALSE]
}
