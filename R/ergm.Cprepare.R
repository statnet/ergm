#  File ergm/R/ergm.Cprepare.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
##########################################################################
# The <ergm.Cprepare> function builds an object called Clist that contains
# all the necessary ingredients to be passed to the C functions
##########################################################################

ergm.Cprepare <- function(nw, m)
{
  n <- network.size(nw)
  dir <- is.directed(nw)
  Clist<-list(n=n, dir=dir)
  bip <- nw$gal$bipartite
  if (is.null(bip)) bip <- 0
  Clist$bipartite <- bip
  Clist$ndyads <- n * (n-1) / (2-dir)
  e<-as.edgelist(nw) # Ensures that for undirected networks, tail<head.
  if(length(e)==0){
    Clist$nedges<-0
    Clist$tails<-NULL
    Clist$heads<-NULL
  }else{
    if(!is.matrix(e)){e <- matrix(e, ncol=2)}
    
    Clist$nedges<-dim(e)[1]
    Clist$tails<-e[,1]
    Clist$heads<-e[,2]
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
                                 } else {
                                   "ergm"
                                 } )
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
ergm.el.lasttoggle <- function(nw){
  n <- network.size(nw)
  b <- if(is.bipartite(nw)) nw %n% "bipartite"
  edge.to.pos <-
    if(is.bipartite(nw))
      function(e) (e[2]-b-1)*b + e[1]
    else if(is.directed(nw))
      function(e) (e[2]-1)*(n - 1) + e[1] - (e[1]>e[2])
    else function(e) (e[2] - 1)*(e[2] - 2)/2 + e[1]

  el <- as.edgelist(nw)
  cbind(el,(nw %n% "lasttoggle")[apply(el,1,edge.to.pos)])
}
