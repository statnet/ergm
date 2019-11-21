#  File R/ergm.Cprepare.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################


#' Internal Functions to Prepare Data for ergm's C Interface
#' 
#' These are internal functions not intended to be called by end
#' users. `ergm_Clist` collates the information in the given object
#' into a form suitable for being passed to the C routines.
#'
#' @param object object to be collated.
#' @param ... additional arguments for methods.
#' @return A list of class `"ergm_Clist"` and possibly a subclass `"ORIGINAL.ergm_Clist"` containing some subset of the following elements: 
#' @keywords internal
#' @export
ergm_Clist <- function(object, ...){
  UseMethod("ergm_Clist")
}

#' @rdname ergm_Clist
#' 
#' @description The \code{ergm.Cprepare} is a legacy function that constructs a combination of `ergm_Clist`s from the given [`network`] and the given [`ergm_model`].
#'
#' @param nw a network or similar object
#' @param m a model object, as returned by \code{\link{ergm_model}}
#' @param verbose logical, whether the design matrix should be printed;
#' default=FALSE
#'
#' @export ergm.Cprepare
ergm.Cprepare <- function(nw, m, response=NULL){
  nw.Clist <- ergm_Clist(nw, response=response)
  m.Clist <- ergm_Clist(m)

  c(nw.Clist, m.Clist)
}


#' @describeIn ergm_Clist
#'
#' Collates a [`network`] object.
#'
#' @template response
#' 
#' @return
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
#' 
#' @export
ergm_Clist.network <- function(object, response=NULL, ...){
  NVL(response) <- NVL(object %ergmlhs% "response")
  e <- na.omit(as.edgelist(object,attrname=response)) # Ensures that for undirected networks, tail<head.
  class(object) <- "network"

  n <- network.size(object)
  dir <- is.directed(object)
  Clist<-list(n=n, dir=dir)
  bip <- object %n% "bipartite"
  if (is.null(bip)) bip <- 0
  Clist$bipartite <- bip
  Clist$ndyads <- network.dyadcount(object)

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

  Clist$lasttoggle <- object %n% "lasttoggle"
  Clist$time <- object %n% "time"

  class(Clist) <- c("network.ergm_Clist", "ergm_Clist")
  Clist
}

#' @noRd
ergm_Clist.pending_update_network <- ergm_Clist.network

#' @describeIn ergm_Clist
#'
#' Collates an [`ergm_model`] object.
#'
#' @return 
#' \item{nterms}{ the number of model terms }
#' \item{nstats}{ the total number of change statistics for all model terms }
#' \item{inputs}{ the concatenated vector of 'input's from each model term as returned by
#' `InitErgmTerm.X` or `InitErgm.X` }
#' \item{fnamestring}{ the concatenated string of model term names }
#' \item{snamestring}{ the concatenated string of package names that contain the C function 'd_fname'; default="ergm" for each fname in fnamestring }
#' @export
ergm_Clist.ergm_model <- function(object, ...){
  mo<-object$terms 
  Clist <- list()

  Clist$nterms<-length(mo)
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
  mo <- object$model.aux$terms
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

  Clist$slots.extra.aux <- unlist(object$slots.extra.aux)
  
  while (substring(Clist$fnamestring, 1, 1)==" ")
    Clist$fnamestring <- substring(Clist$fnamestring, 2)
  while (substring(Clist$snamestring, 1, 1)==" ")
    Clist$snamestring <- substring(Clist$snamestring, 2)

  
  
  # We don't care about diagnostics for terms that are not being
  # estimated.
  Clist$diagnosable <- ! object$etamap$offsetmap
  names(Clist$diagnosable) <- object$coef.names
    
  class(Clist) <- c("ergm_model.ergm_Clist", "ergm_Clist")
  Clist
}

#' @describeIn to_ergm_Cdouble
#'
#' Method for [`network`] objects.
#'
#' @param attrname name of an edge attribute.
#' 
#' @export
to_ergm_Cdouble.network <- function(x, attrname=NULL, ...){
  xm <- as.edgelist(x, attrname=attrname)
  c(nrow(xm),c(na.omit(xm)))
}

#' @noRd
to_ergm_Cdouble.pending_update_network <- to_ergm_Cdouble.network


#' @describeIn to_ergm_Cdouble
#'
#' Method for [`matrix`] objects, assumed to be edgelists.
#'
#' @param prototype A network whose relevant attributes (size,
#'   directedness, bipartitedness, and presence of loops) are imposed
#'   on the output edgelist if \code{x} is already an edgelist. (For
#'   example, if the prototype is undirected, `to_ergm_Cdouble`
#'   will ensure that \eqn{t < h}.)
#' @keywords internal
#' @export
to_ergm_Cdouble.matrix <- function(x, prototype=NULL, ...){
  x <- if(!is.null(prototype)) as.edgelist(x, n=network.size(prototype), directed=is.directed(prototype),
                                           bipartite=if(is.bipartite(prototype)) prototype%n%"bipartite" else 0,
                                           loops=has.loops(prototype))
       else x[order(x[,1],x[,2]),,drop=FALSE]
  c(nrow(x),c(na.omit(x)))
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
#' Though this is an API intended for internal use, some functions,
#' like \code{\link[tergm]{stergm}} (for EGMME),
#' \code{\link[tergm:simulate.stergm]{simulate}}, and
#' \code{\link[=summary.formula]{summary}} can be passed networks with
#' this information using the following \code{\link{network}} (i.e.,
#' \code{\link{\%n\%}}) attributes:
#' \describe{
#'
#' \item{`"time"`}{the
#' time stamp associated with the network}
#'
#' \item{`"lasttoggle"`}{an integer vector analogous to the one returned by [to_ergm_Cdouble.network()] with an attribute: number of elements, a list of tails, a list of heads, and time the edge was last toggled.}
#' }
#'
#' On the C side, it is represented by a hash table.
#' 
#' Again, this API is subject to change without notice.
#'
#' @aliases lasttoggle last.toggle last-toggle
#' @keywords internal
#' @name lasttoggle
NULL
