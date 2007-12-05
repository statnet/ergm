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
### This section commented out by DH on Oct. 13, 2007.  Do we really need it?
### If so, note that adding it back in will introduce a bug in 
### the MPLE routine, which already knows how to deal with bipartite
### networks without being passed a "design" matrix.  This Clist.miss only 
### confuses the MPLE_wrapper routine.
#### MSH on Oct 13 notes that this function (ergm.design) was added in rev 47
#### to deal with missing data for bipartite networks.
#### It is clearly no longer needed.
#    if(is.bipartite(nw)){                       
#     nb1 <- get.network.attribute(nw,"bipartite")
#     nb2 <- network.size(nw) - nb1
#     temp <- matrix(0,ncol=nb1,nrow=nb1)
#     base <- cbind(as.vector(col(temp)), as.vector(row(temp)))
#     base <- base[base[, 2] > base[, 1], ]
#     if(Clist.miss$nedges==0){
#      Clist.miss$heads <- NULL
#      Clist.miss$tails <- NULL
#     }
#     Clist.miss$heads <- c(Clist.miss$heads,base[,2]) 
#     Clist.miss$tails <- c(Clist.miss$tails,base[,1]) 
#     temp <- matrix(0,ncol=nb2,nrow=nb2)
#     base <- cbind(as.vector(col(temp)), as.vector(row(temp)))
#     base <- base[base[, 2] > base[, 1], ]
#     Clist.miss$heads <- c(Clist.miss$heads,base[,2]+nb1) 
#     Clist.miss$tails <- c(Clist.miss$tails,base[,1]+nb1) 
#     Clist.miss$nedges<-Clist.miss$nedges+(nb1*(nb1-1)+nb2*(nb2-1))/2
#    }
   }
   Clist.miss
}
