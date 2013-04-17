#  File R/ergm.Cprepare.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
##########################################################################
# The <ergm.Cprepare> function builds an object called Clist that contains
# all the necessary ingredients to be passed to the C functions
#
# --PARAMETERS--
#   nw:  a network object
#   m :  a model object, as returned by <ergm.getmodel>
#
# --RETURNED--
#   Clist:  a list of parameters used by several of the fitting routines
#           containing
#            n           :  the size of the network
#            dir         :  whether the network is directed (T or F)
#            bipartite   :  whether the network is bipartite (T or F)
#            ndyads      :  the number of dyads in the network
#            nedges      :  the number of edges in this network
#            tails       :  the vector of tail nodes; tail nodes are
#                               the 1st column of the implicit edgelist,
#                               so either the lower-numbered nodes in an
#                               undirected graph, or the out nodes of a
#                               directed graph, or the b1 nodes of a bi-
#                               partite graph
#            heads           :  the vector of head nodes; head nodes are
#                               the 2nd column of the implicit edgelist,
#                               so either the higher-numbered nodes in an
#                               undirected graph, or the in nodes of a
#                               directed graph, or the b2 nodes of a bi-
#                               partite graph
#            nterms      :  the number of model terms
#            nstats      :  the total number of change statistics
#                           for all model terms
#            inputs      :  the concatenated vector of 'input's from each
#                           model term as returned by <InitErgmTerm.X> or
#                           <InitErgm.X>
#            fnamestring :  the concatenated string of model term names
#            snamestring :  the concatenated string of package names that
#                           contain the C function 'd_fname'; default="ergm"
#                           for each fname in fnamestring
#
##########################################################################

ergm.Cprepare <- function(nw, m, response=NULL)
{
  n <- network.size(nw)
  dir <- is.directed(nw)
  Clist<-list(n=n, dir=dir)
  bip <- nw$gal$bipartite
  if (is.null(bip)) bip <- 0
  Clist$bipartite <- bip
  Clist$ndyads <- n * (n-1) / (2-dir)
  e<-as.edgelist(nw,attrname=response) # Ensures that for undirected networks, tail<head.
  if(length(e)==0){
    Clist$nedges<-0
    Clist$tails<-NULL
    Clist$heads<-NULL
    ## Make sure weights is not NULL if response!=NULL, even if it's
    ## empty, since it's used to decide whether MCMC or WtMCMC is
    ## called.
    if(!is.null(response)) Clist$weights<-numeric(0)
  }else{
    if(!is.matrix(e)){e <- matrix(e, ncol=2+!is.null(response))}

    ## Delete 0 edges.
    if(!is.null(response)) e<-e[e[,3]!=0,,drop=FALSE]
    
    Clist$nedges<-dim(e)[1]
    Clist$tails<-e[,1]
    Clist$heads<-e[,2]
    if(!is.null(response)) Clist$weights<-e[,3]
  }

  Clist$lasttoggle <- nw %n% "lasttoggle"
  Clist$time <- nw %n% "time"
  
  mo<-m$terms 
  
  Clist$nterms<-length(mo)
  Clist$nstats<-0
  Clist$fnamestring<-""
  Clist$snamestring<-""
  Clist$inputs<-numeric(0)
  if (Clist$nterms>0) {
    for(i in 1:Clist$nterms) {
      term_i <- mo[[i]]
      Clist$fnamestring <- paste(Clist$fnamestring, term_i$name)
      # This lets "pkgname" play the same role as "soname":
      Clist$snamestring <- paste(Clist$snamestring, 
                                 if (!is.null(term_i$soname)) {
                                   term_i$soname
                                 } else if (!is.null(term_i$pkgname)) {
                                   term_i$pkgname
                                 } else stop("ERGM term specifying C function `", term_i$name,"' is missing C library or package name.") )
      Clist$inputs <- c(Clist$inputs, term_i$inputs)
      Clist$nstats <- Clist$nstats + term_i$inputs[2]
    }
  }
  while (substring(Clist$fnamestring, 1, 1)==" ")
    Clist$fnamestring <- substring(Clist$fnamestring, 2)
  while (substring(Clist$snamestring, 1, 1)==" ")
    Clist$snamestring <- substring(Clist$snamestring, 2)

  # We don't care about diagnostics for terms that are not being
  # estimated.
  Clist$diagnosable <- ! m$etamap$offsetmap
  names(Clist$diagnosable) <- m$coef.names[!m$etamap$offsetmap]
    
  Clist
}


## Construct and serialize a very simple static edgelist, with the
## vertex having the lesser index the tail and sorted by tails, then
## by heads.
ergm.Cprepare.el<-function(x, attrname=NULL, directed=if(is.network(x)) is.directed(x) else stop("Directedness argument is mandatory for edgelist input.")){
  xm <- if(is.network(x)) as.edgelist(x, attrname=attrname) else x
  
  if(nrow(xm)){
    # Sort.
    xm <- xm[order(xm[,1],xm[,2]),,drop=FALSE]
  }

  c(nrow(xm),c(xm))
}

# Note: this converter must be kept in sync with whatever edgetree.c does.
mk.edge.to.pos.lasttoggle.f <- function(nw){
  if(is.bipartite(nw)){
    b <- if(is.bipartite(nw)) nw %n% "bipartite"
    function(e) (e[2] - b - 1)*b + e[1]
  }else{
    n <- network.size(nw)
    if(is.directed(nw))
      function(e) (e[2] - 1)*(n - 1) + e[1] - (e[1] > e[2])
    else
      function(e) (e[2] - 1)*(e[2] - 2)/2 + e[1]
  }
}

ergm.el.lasttoggle <- function(nw){
  edge.to.pos <- mk.edge.to.pos.lasttoggle.f(nw)
  el <- as.edgelist(nw)
  cbind(el,(nw %n% "lasttoggle")[apply(el,1,edge.to.pos)])
}

to.matrix.lasttoggle <- function(nw){
  n <- network.size(nw)
  b <- if(is.bipartite(nw)) nw %n% "bipartite"
  
  if(is.bipartite(nw)) m <- matrix(nw %n% "lasttoggle", b, n-b, byrow=FALSE)
  else{
    m <- matrix(0,n,n)
    if(is.directed(nw))
      m[as.logical(1-diag(1,nrow=n))] <- nw %n% "lasttoggle"
    else{      
      m[upper.tri(m)] <- nw %n% "lasttoggle"
      m <- m + t(m)
    }
  }
  m
}

to.lasttoggle.matrix <- function(m, directed=TRUE, bipartite=FALSE){
  if(bipartite) c(m)
  else if(directed) c(m[as.logical(1-diag(1,nrow=nrow(m)))])
  else c(m[upper.tri(m)])
}
