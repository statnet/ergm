ergm.design <- function(nw, model, verbose=FALSE){
   notobserved <- is.na(nw)
   if(network.edgecount(notobserved)==0){
    Clist.miss <- list(heads=NULL, tails=NULL, nedges=0, dir=is.directed(nw))
   }else{
    Clist.miss <- ergm.Cprepare(notobserved, model)
    if(verbose){
      cat("Design matrix:\n")
      summary(notobserved)
    }
   }
   Clist.miss
}
