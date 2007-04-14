ergm.design <- function(nw, model, initialfit=FALSE, verbose=FALSE){
   notobserved <- get.network.attribute(nw,"design")
   if(is.null(notobserved)){
    mClist <- list(heads=0, tails=0, nedges=0, dir=is.directed(nw))
   }else{
    mClist <- ergm.Cprepare(notobserved, model)
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
     if(mClist$nedges==0){
      mClist$heads <- NULL
      mClist$tails <- NULL
     }
     mClist$heads <- c(mClist$heads,base[,2]) 
     mClist$tails <- c(mClist$tails,base[,1]) 
     temp <- matrix(0,ncol=nevents,nrow=nevents)
     base <- cbind(as.vector(col(temp)), as.vector(row(temp)))
     base <- base[base[, 2] > base[, 1], ]
     mClist$heads <- c(mClist$heads,base[,2]+nactors) 
     mClist$tails <- c(mClist$tails,base[,1]+nactors) 
     mClist$nedges<-mClist$nedges+(nactors*(nactors-1)+nevents*(nevents-1))/2
    }
   }
   mClist
}
