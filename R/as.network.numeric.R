#  File ergm/R/as.network.numeric.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
###########################################################################
# The <as.network.numeric> function creates and returns a bernouli
# network.
##########################################################################
as.network.numeric<-function(x,
    directed = TRUE,
    hyper = FALSE, loops = FALSE, multiple = FALSE, bipartite = FALSE,
    ignore.eval = TRUE, names.eval = NULL,
    edge.check = FALSE,
    density=NULL, init=NULL, numedges=NULL, ...){
  #returns a bernouli network.
  if(bipartite){
   nb2 <- x
   nb1 <- bipartite
   directed <- FALSE
  }else{                                                        
   nb2 <- x
   nb1 <- x
  }
  if(directed)
    ndyads <- nb1*(nb1-1)
  else if(bipartite)
    ndyads <- nb1*nb2
  else
    ndyads <- nb1*(nb1-1)/2
  
  if(missing(density)){
    if(missing(init)){
      #     So the expected number of ties is the same as
      #     the number of nodes
      density <- nb1/ndyads
    }else{
      density <- exp(init)/(1+exp(init))
    }
  }
  nw.mat <- matrix(0,nrow=nb1,ncol=nb2)
  dimnames(nw.mat) <- list(1:nb1,1:nb2)
  if(is.null(numedges)){
    nwij <- runif(ndyads)<density
  }else{
    nwij <- rep(0,ndyads)
    nwij[sample(1:ndyads,size=numedges,replace=FALSE)] <- 1
  }
  if(directed)
    nw.mat[row(nw.mat)!=col(nw.mat)] <- nwij
  else if(bipartite)
    nw.mat[,] <- nwij
  else{
    nw.mat[row(nw.mat) < col(nw.mat)] <- nwij
    nw.mat <- nw.mat + t(nw.mat)
  }
  #Return the result
  network(nw.mat,directed=directed, bipartite=bipartite>0)
}
