as.edgelist <- function(nw, attrname = NULL, as.sna.edgelist = FALSE,...){
  el <- as.matrix.network.edgelist(nw, attrname=attrname, as.sna.edgelist=as.sna.edgelist,...)
  if(!is.directed(nw)) el[,1:2] <- cbind(pmin(el[,1],el[,2]),pmax(el[,1],el[,2]))
  el[order(el[,1],el[,2]),,drop=FALSE]
}
