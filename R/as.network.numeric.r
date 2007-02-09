# Added by MSH 4/9/06 to allow the simple formation of networks
as.network.numeric<-function(x,
    directed = TRUE,
    hyper = FALSE, loops = FALSE, multiple = FALSE, bipartite = FALSE,
    ignore.eval = TRUE, names.eval = NULL,
    edge.check = FALSE,
    density=NULL, theta0=NULL, numedges=NULL, ...){
  #returns a bernouli network.
  if(bipartite){
   nevents <- x
   nactors <- bipartite
   directed <- FALSE
  }else{
   nevents <- x
   nactors <- x
  }
  if(directed)
    ndyads <- nactors*(nactors-1)
  else if(bipartite)
    ndyads <- nactors*nevents
  else
    ndyads <- nactors*(nactors-1)/2
  
  if(missing(density)){
    if(missing(theta0)){
      #     So the expected number of ties is the same as
      #     the number of nodes
      density <- nactors/ndyads
    }else{
      density <- exp(theta0)/(1+exp(theta0))
    }
  }
  nw.mat <- matrix(0,nrow=nactors,ncol=nevents)
  dimnames(nw.mat) <- list(1:nactors,1:nevents)
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
