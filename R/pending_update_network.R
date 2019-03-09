#  File R/pending_update_network.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################
#' A (Relatively) Lightweight Read-Only Representation of Network Objects
#'
#' `pending_update_network` is an semi-internal class for passing
#' around results of MCMC sampling, particularly when the result is
#' used to start another MCMC sampler. It is deliberately loosely
#' specified, and its structure and even name are subject to change.
#'
#' @param nw a [`network`] object.
#' @param update a character string, a list with elements named
#'   `"newnwtails"`, `"newnwheads"`, and (optionally)
#'   `"newnwweights"`, or `"newedgelist"`, or `NULL`. See Details.
#' @param response,attrname Name of edge attribute to get or set. If
#'   `NULL`, binary network is assumed.
#'
#' @return
#'
#' At this time, a `pending_update_network` object is (subject to
#' change) a [`network`] object with all edges removed and with a
#' network attribute `".update"` containing a two or three column
#' matrix with a (possibly valued) edge list. The third column name
#' also indicates the edge attribute represented.
#'
#' If `update` uses `"newedgelist"`, it is copied directly; othewirse,
#' `"newnwtails"`, `"newnwheads"`, and (optionally) `"newnwweights"`
#' are converted to an edgelist.
#'
#' Its class is set (not subclassed!) to `pending_update_network`, in
#' order to prevent [`network`] accessors and modifiers from affecting
#' it.
#'
#' @details The `update` argument controls the contents of the object
#'   depending on its mode: \describe{
#'
#' \item{`NULL`}{Use `nw`'s own edgelist, and don't set weights.}
#'
#' \item{`character`}{Use `nw`'s own edgelist, and don't set weights
#' for attribute specified by `update`.}
#'
#' \item{`a list`}{Use the elements of the list, dropping others.}
#'
#' }
#'
#' @keywords internal
#' @export
pending_update_network <- function(nw, update=NULL, response=if(is.character(update)) update){
  if(!is.network(nw) && !is.pending_update_network(nw)) stop("nw must be a network or a pending_update_network object.")

  el <-
    # Note that the following will get the edgelist from a pending_update_network as well.
    if(is.null(update)||is.character(update)) as.edgelist(nw, attrname=response)
    else .extract_z_edgelist(update, response=response)

  if(!is.null(response)&&ncol(el)==3) colnames(el)[3]<-response

  if(is(nw, "pending_update_network")) class(nw) <- "network"
  nw[,] <- FALSE
  nw%n%".update" <- el
  
  class(nw) <- "pending_update_network"
  nw
}

#' @rdname pending_update_network
#' @export
is.pending_update_network <- function(x){
  is(x, "pending_update_network")
}

#' @rdname pending_update_network
#' @export
as.edgelist.pending_update_network <- function(x,...){
  class(x) <- "network"
  e <- x%n%".update"
  if(length(e)!=0) e <- e[order(e[,1],e[,2]),,drop=FALSE]
  attr(e, "n") <- network.size(x)
  attr(e, "vnames") <- x%v%"vertex.names"
  attr(e, "directed") <- is.directed(x)
  attr(e, "bipartite") <- x%n%"bipartite"
  attr(e, "loops") <- has.loops(x)
  e
}

#' @rdname pending_update_network
#' @export
as.matrix.pending_update_network <- function(x,matrix.type=NULL,attrname=NULL,...){
  matrix.type<-match.arg(matrix.type,c("adjacency","incidence","edgelist"))
  if(matrix.type!="edgelist") stop("pending_update_network can only be converted to an edgelist at this time")

  x <- as.network(x,populate=FALSE)

  if(!is.null(attrname) && colnames(x%n%".update")[3]!=attrname) stop("Edge attribute ",sQuote(attrname)," is not stored in the pending_update_network object.")
  x %n% ".update"
}

.extract_z_edgelist <- function(z, response=NULL){
  # if z has a newedgelist attached, use it
  if("newedgelist" %in% names(z)){
    newedgelist<-z$newedgelist[,1:2,drop=FALSE]
    newnwweights<- if(ncol(z$newedgelist)==3) z$newedgelist[,3]
  }else{
    # expect that z will have seperate lists of heads and tails
    nedges<-z$newnwtails[1]
    # *** don't forget - edgelists are cbind(tails, heads) now
    newedgelist <- cbind(z$newnwtails[seq_len(nedges)+1],z$newnwheads[seq_len(nedges)+1])
    newnwweights <- z$newnwweights[seq_len(nedges)+1]
  }
  out <- cbind(newedgelist, newnwweights)
  if(!is.null(response)&&ncol(out)==3) colnames(out)[3]<-response
  out
}

#' @rdname pending_update_network
#' @export
as.network.pending_update_network <- function(x, ..., populate=TRUE){
  if(!populate){
    class(x) <- c("pending_update_network","network")
    return(x)
  }

  class(x) <- c("network")

  z <- x%n%".update"
  delete.network.attribute(x, ".update")

  update(x,z,matrix.type="edgelist", attrname=if(ncol(z)==3) colnames(z)[3])
}

#' @describeIn pending_update_network Note that this method fails when
#'   `na.omit=FALSE`, since missing edges are not stored.
#' @param na.omit Whether missing edges should be counted. Note that
#'   missing edge information is not stored.
#' @export
network.edgecount.pending_update_network <- function(x, na.omit=TRUE,...){
  if(!na.omit) stop("pending_update_network cannot store missing edges.")
  nrow(as.network(x,populate=FALSE)%n%".update")
}

#' @describeIn pending_update_network Note that this method fails with
#'   its default argument, since missing edges are not stored.
#' @export
network.dyadcount.pending_update_network <- function(x, na.omit=TRUE,...){
  if(na.omit) stop("pending_update_network cannot store missing edges.")
  class(x) <- "network"
  network.dyadcount(x)
}

#' @rdname pending_update_network
#' @export
network.size.pending_update_network <- function(x,...){
  NextMethod()
}

#' @describeIn pending_update_network A stub that produces an error.
#' @export
network.naedgecount.pending_update_network <- function(x,...){
  stop("pending_update_network cannot store missing edges.")
}
