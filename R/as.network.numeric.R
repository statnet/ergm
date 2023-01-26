#  File R/as.network.numeric.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2023 Statnet Commons
################################################################################

#' Create a Simple Random network of a Given Size
#' 
#' \code{\link{as.network.numeric}} creates a random Bernoulli network of the
#' given size as an object of class \code{\link[network]{network}}.
#' 
#' The network will not have vertex, edge or network attributes.  These
#' can be added with operators such as \code{\%v\%}, \code{\%n\%}, \code{\%e\%}.
#' 
#' @param x count; the number of nodes in the network
#' @param directed logical; should edges be interpreted as directed?
#' @param hyper logical; are hyperedges allowed? Currently ignored.
#' @param loops logical; should loops be allowed? Currently ignored.
#' @param multiple logical; are multiplex edges allowed? Currently ignored.
#' @param bipartite count; should the network be interpreted as bipartite? If
#' present (i.e., non-NULL) it is the count of the number of actors in the
#' bipartite network. In this case, the number of nodes is equal to the number
#' of actors plus the number of events (with all actors preceding all events).
#' The edges are then interpreted as nondirected.
#' @param ignore.eval logical; ignore edge values? Currently ignored.
#' @param names.eval optionally, the name of the attribute in which edge values
#' should be stored. Currently ignored.
#' @param edge.check logical; perform consistency checks on new edges?
#' @param density numeric; the probability of a tie for Bernoulli networks. If
#' neither density nor \code{init} is given, it defaults to the number of nodes
#' divided by the number of dyads (so the expected number of ties is the same
#' as the number of nodes.)
#' @param init numeric; the log-odds of a tie for Bernoulli networks.  It is
#' only used if density is not specified.
#' @param numedges count; if present, sample the Bernoulli network conditional
#' on this number of edges (rather than independently with the specified
#' probability).
#' @param ... additional arguments
#' @return An object of class \code{\link[network]{network}}
#' @seealso \code{\link[network]{network}}
#' @references Butts, C.T.  2002.  ``Memory Structures for Relational Data in
#' R: Classes and Interfaces'' Working Paper.
#' @keywords classes graphs
#' @examples
#' 
#' #Draw a random directed network with 25 nodes
#' g<-network(25)
#' #Draw a random undirected network with density 0.1
#' g<-network(25, directed=FALSE, density=0.1)
#' #Draw a random bipartite network with 4 actors and 6 events and density 0.1
#' g<-network(10, bipartite=4, density=0.1)
#' @importFrom network as.network
#' @export
as.network.numeric<-function(x,
    directed = TRUE,
    hyper = FALSE, loops = FALSE, multiple = FALSE, bipartite = FALSE,
    ignore.eval = TRUE, names.eval = NULL,
    edge.check = FALSE,
    density=NULL, init=NULL, numedges=NULL, ...){
  #returns a bernouli network.
  if(bipartite){
   nb2 <- x - bipartite
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
