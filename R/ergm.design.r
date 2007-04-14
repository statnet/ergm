ergm.design <- function(nw, model, initialfit=FALSE, verbose=FALSE){
   notobserved <- get.network.attribute(nw,"design")
   if(is.null(notobserved)){
    Clist.miss <- list(heads=0, tails=0, nedges=0, dir=is.directed(nw))
   }else{
    Clist.miss <- ergm.Cprepare(notobserved, model)
    if(verbose){
     cat("Design matrix:\n")
     summary(notobserved)
    }
   }
   if(initialfit){
    if(is.bipartite(nw)){
     nactors <- get.network.attribute(nw,"bipartite")
     nevents <- network.size(nw) - nactors
     temp <- matrix(0,ncol=nactors,nrow=nactors)
     base <- cbind(as.vector(col(temp)), as.vector(row(temp)))
     base <- base[base[, 2] > base[, 1], ]
     if(Clist.miss$nedges==0){
      Clist.miss$heads <- NULL
      Clist.miss$tails <- NULL
     }
     Clist.miss$heads <- c(Clist.miss$heads,base[,2]) 
     Clist.miss$tails <- c(Clist.miss$tails,base[,1]) 
     temp <- matrix(0,ncol=nevents,nrow=nevents)
     base <- cbind(as.vector(col(temp)), as.vector(row(temp)))
     base <- base[base[, 2] > base[, 1], ]
     Clist.miss$heads <- c(Clist.miss$heads,base[,2]+nactors) 
     Clist.miss$tails <- c(Clist.miss$tails,base[,1]+nactors) 
     Clist.miss$nedges<-Clist.miss$nedges+(nactors*(nactors-1)+nevents*(nevents-1))/2
    }
   }
   Clist.miss
}
