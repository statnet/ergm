#  File R/as.network.numeric.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
###########################################################################
# The <as.network.numeric> function creates and returns a bernouli
# network.
#
# --PARAMETERS--
#   x        : for a non-bipartite network, the number of nodes;
#              for a bipartite network, the number of events.
#               (the number of actors is specied via the bipartite param)
#   directed : whether the network is to be directed; default=TRUE
#   bipartite: the count of actors if the network should be bipartite; 0
#              if 'x' is not bipartite; default=FALSE
#   density  : the probability of a tie; default=the number of nodes divided
#              by the number of possible dyad IF init isn't provided, NULL
#              otherwise
#   init   : the log-odds of a tie, this parameter is ignored if density
#              is given; default=the number of nodes divided by the number of
#              possible dyad IF density isn't provided, NULL otherwise
#   numedges : the number of edges that the returned network must have;
#              default=NULL, in which case numedges will result from
#              the random process
#
#
# --IGNORED PARAMETERS--
#   hyper       : whether the network should allow hyper edges; default=FALSE
#   loops       : whether the network should allow loops; default=FALSE
#   multiple    : whether the network should allow multiplex edges;
#                 default=FALSE
#   ignore.eval : whether edge values should be ignored; default=FALSE
#                 default=FALSE
#   names.eval  : the attribute in which edge values are to be stored;
#                 default=NULL
#   edge.check  : whether a consistency check should be performed;
#                 default=FALSE
#   ...         : additional parameters
#
#
# --RETURNED--
#   a random bernoulli network with the specified size and desired
#   probabilistic qualities
#
# author: MSH
#
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
