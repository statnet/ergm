#  File R/ergm.Cprepare.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################

#' Internal Function to Prepare Data for ergm's C Interface
#' 
#' These are internal functions not intended to be called by end users.  The
#' \code{ergm.Cprepare} function builds an object called `Clist` that contains
#' all the necessary ingredients to be passed to the C functions, other
#' functions create edgelists and handle missing edge data. These low-level functions are used by other ergm-related packages, but
#' should never need to be called directly by the user.
#'
#' @param nw,x a network or similar object
#' @param m a model object, as returned by \code{\link{ergm.getmodel}}
#' @template response
#' @param verbose logical, whether the design matrix should be printed;
#' default=FALSE
#' @return \code{ergm.Cprepare} returns `Clist`: a list of parameters used by
#' several of the fitting routines containing
#' \item{n}{ the size of the network }
#' \item{dir}{ whether the network is directed (T or F) }
#' \item{bipartite}{ whether the network is bipartite (T or F) }
#' \item{ndyads}{ the number of dyads in the network }
#' \item{nedges}{ the number of edges in this network }
#' \item{tails}{ the vector of tail nodes; tail nodes are the 1st
#' column of the implicit edgelist, so either the lower-numbered nodes in an
#' undirected graph, or the out nodes of a directed graph, or the b1 nodes of a
#' bipartite graph }
#' \item{heads}{ the vector of head nodes; head nodes are the
#' 2nd column of the implicit edgelist, so either the higher-numbered nodes in
#' an undirected graph, or the in nodes of a directed graph, or the b2 nodes of
#' a bipartite graph }
#' \item{nterms}{ the number of model terms }
#' \item{nstats}{ the total number of change statistics for all model terms }
#' \item{inputs}{ the concatenated vector of 'input's from each model term as returned by
#' `InitErgmTerm.X` or `InitErgm.X` }
#' \item{fnamestring}{ the concatenated string of model term names }
#' \item{snamestring}{ the concatenated string of package names that contain the C function 'd_fname'; default="ergm" for each fname in fnamestring }
#' 
#' @export ergm.Cprepare
ergm.Cprepare <- function(nw, m, response=NULL)
{
  e<-as.edgelist(nw,attrname=response) # Ensures that for undirected networks, tail<head.
  class(nw) <- "network"

  n <- network.size(nw)
  dir <- is.directed(nw)
  Clist<-list(n=n, dir=dir)
  bip <- nw %n% "bipartite"
  if (is.null(bip)) bip <- 0
  Clist$bipartite <- bip
  Clist$ndyads <- network.dyadcount(nw)

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
  
  Clist$nstats<-0
  Clist$fnamestring<-""
  Clist$snamestring<-""
  Clist$inputs<-numeric(0)
  
  Clist$nterms<-length(mo)
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

  # Attach the auxiliaries
  mo <- m$model.aux$terms
  anterms <- length(mo)
  Clist$nterms <- Clist$nterms + anterms 
  if (anterms>0) {
    for(i in 1:anterms) {
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
      # Auxiliaries do not produce stats.
    }
  }

  Clist$slots.extra.aux <- unlist(m$slots.extra.aux)
  
  while (substring(Clist$fnamestring, 1, 1)==" ")
    Clist$fnamestring <- substring(Clist$fnamestring, 2)
  while (substring(Clist$snamestring, 1, 1)==" ")
    Clist$snamestring <- substring(Clist$snamestring, 2)

  
  
  # We don't care about diagnostics for terms that are not being
  # estimated.
  Clist$diagnosable <- ! m$etamap$offsetmap
  names(Clist$diagnosable) <- m$coef.names
    
  Clist
}



#' @rdname ergm.Cprepare
#' @description `ergm.Cprepare.el` constructs and serializes a very simple static
#'   edgelist, with the vertex having the lesser index the tail and
#'   sorted by tails, then by heads.
#' @param prototype A network whose relevant attributes (size,
#'   directedness, bipartitedness, and presence of loops) are imposed
#'   on the output edgelist if \code{x} is already an edgelist. (For
#'   example, if the prototype is undirected, \code{ergm.Cprepare.el}
#'   will ensure that \eqn{t < h}.)
#' @param attrname name of an edge attribute.
#' @export ergm.Cprepare.el
ergm.Cprepare.el<-function(x, attrname=NULL, prototype=NULL){
  xm <- if(is.network(x) || is(x, "pending_update_network")) as.edgelist(x, attrname=attrname)
        else if(is(x, "rlebdm")) as.edgelist(x, prototype=prototype)
        else if(!is.null(prototype)) as.edgelist.matrix(x, n=network.size(prototype), directed=is.directed(prototype),
                                                        bipartite=if(is.bipartite(prototype)) prototype%n%"bipartite" else 0,
                                                        loops=has.loops(prototype))
        else x[order(x[,1],x[,2]),,drop=FALSE]
                                                        
  c(nrow(xm),c(xm))
}

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

#' Storing last toggle information in a network
#' 
#' An informal extension to \code{\link{network}} objects allowing
#' some limited temporal information to be stored.
#' WARNING: THIS DOCUMENTATION IS PROVIDED AS A COURTESY, AND THE API
#' DESCRIBED IS SUBJECT TO CHANGE WITHOUT NOTICE, DOWN TO COMPLETE
#' REMOVAL. NOT ALL FUNCTIONS THAT COULD SUPPORT IT DO. USE AT YOUR
#' OWN RISK.
#' 
#' While \code{\link[networkDynamic]{networkDynamic}} provides a flexible,
#' consistent method for storing dynamic networks, the \code{C} routines of
#' \code{\link[=ergm-package]{ergm}} and
#' \code{\link[tergm:tergm-package]{tergm}} required a simpler and more
#' lightweight representation.
#' 
#' This representation consisted of a single integer representing the
#' time stamp and an integer vector of length to
#' \code{\link{network.dyadcount}(nw)} --- the number of potential
#' ties in the network, giving the last time point during which each
#' of the dyads in the network had changed.
#' 
#' Though this is an API intended for internal use, some functions,
#' like \code{\link[tergm]{stergm}} (for EGMME),
#' \code{\link[tergm:simulate.stergm]{simulate}}, and
#' \code{\link[=summary.formula]{summary}} can be passed networks with
#' this information using the following \code{\link{network}} (i.e.,
#' \code{\link{\%n\%}}) attributes: \describe{ \item{list("time")}{the
#' time stamp associated with the network} \item{list("lasttoggle")}{a
#' vector of length \code{\link{network.dyadcount}(nw)}, giving the
#' last change time associated with each dyad. See the source code of
#' \code{\link[=ergm-package]{ergm}} internal functions
#' \code{to.matrix.lasttoggle}, \code{ergm.el.lasttoggle}, and
#' \code{to.lasttoggle.matrix} for how they are serialized.} }
#' 
#' For technical reasons, the \code{\link[tergm:tergm-package]{tergm}}
#' routines treat the \code{lasttoggle} time points as shifted by
#' \eqn{-1}.
#' 
#' Again, this API is subject to change without notice.
#'
#' @aliases lasttoggle last.toggle last-toggle
#' @name lasttoggle
NULL

#' @describeIn lasttoggle Returns a 3-column matrix whose first two
#'   columns are tails and heads of extant edges and whose third
#'   column are the creation times for those edges.
#' @param nw the network, otpionally with a `"lasttoggle"` network
#'   attribute.
#' @export
ergm.el.lasttoggle <- function(nw){
  edge.to.pos <- mk.edge.to.pos.lasttoggle.f(nw)
  el <- as.edgelist(nw)
  cbind(el,NVL((nw %n% "lasttoggle"),0)[apply(el,1,edge.to.pos)]) # change to 0 if null
}

#' @describeIn lasttoggle Returns a numeric sociomatrix whose values
#'   are last toggle times for the corresponding dyads.
#' @export
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

#' @describeIn lasttoggle Serializes a matrix of last toggle times
#'   into the form used by C code.
#' @param m a sociomatrix of appropriate dimension (rectangular for
#'   bipartite networks).
#' @param directed,bipartite whether the matrix represents a directed
#'   and/or a bipartite networks.
#' @export
to.lasttoggle.matrix <- function(m, directed=TRUE, bipartite=FALSE){
  if(bipartite) c(m)
  else if(directed) c(m[as.logical(1-diag(1,nrow=nrow(m)))])
  else c(m[upper.tri(m)])
}
