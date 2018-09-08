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
#' matrix with a (possibly valued) edge list.
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
pending_update_network <- function(nw, update=NULL){
  if(is(nw, "pending_update_network")) class(nw) <- "network"
  if(!is.network(nw)) stop("nw must be a network or a pending_update_network object.")
  
  el <-
    if(is.null(update)||is.character(update)) as.edgelist(nw, attrname=update)
    else .extract_z_edgelist(update, response=TRUE)

  nw <- empty_network(nw)
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
as.edgelist.pending_update_network <- function(x,attrname=NULL,...){
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

.extract_z_edgelist <- function(z, response=NULL){
  # if z has a newedgelist attached, use it
  if("newedgelist" %in% names(z)){
    newedgelist<-z$newedgelist[,1:2,drop=FALSE]
    newnwweights<- if(!is.null(response) && ncol(z$newedgelist)==3) z$newedgelist[,3]
  }else{
    # expect that z will have seperate lists of heads and tails
    nedges<-z$newnwtails[1]
    # *** don't forget - edgelists are cbind(tails, heads) now
    newedgelist <- cbind(z$newnwtails[seq_len(nedges)+1],z$newnwheads[seq_len(nedges)+1])
    newnwweights <- z$newnwweights[seq_len(nedges)+1]
  }
  cbind(newedgelist, newnwweights)
}

#' @rdname pending_update_network
#' @export
as.network.pending_update_network <- function(x, response=NULL, ...){
  class(x) <- "network"
  z <- x%n%".update"
  delete.network.attribute(x, ".update")

  update(x,z,matrix.type="edgelist", attrname=response)
}
