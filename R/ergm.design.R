ergm.design <- function(nw, model, verbose=FALSE){
  if(network.naedgecount(nw)==0){
    Clist.miss <- list(heads=NULL, tails=NULL, nedges=0, dir=is.directed(nw))
  }else{
    Clist.miss <- ergm.Cprepare(is.na(nw), model)
    if(verbose){
      cat("Design matrix:\n")
      print(summary(is.na(nw)))
    }
  }
  Clist.miss
}
