#  File R/ergm_state.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################
#' A Representation of ERGM state
#' 
#' `ergm_state` is an semi-internal class for passing
#' around results of MCMC sampling, particularly when the result is
#' used to start another MCMC sampler. It is deliberately loosely
#' specified, and its structure and even name are subject to change.
#'
#' @param nw a [`network`] object.
#' @param model an [`ergm_model`] object.
#' @param response a character vector representing the response attribute. If
#'   `NULL`, binary network is assumed.
#'
#' @return 
#' At this time, an `ergm_state` object is (subject to
#' change) a list containing the following elements:
#' \describe{
#' 
#' \item{el}{a [`tibble`] [`edgelist`] representing the edge state of the network}
#' 
#' \item{nw0}{a [`network`] object with all edges removed.}
#'
#' \item{stats}{a numeric vector of network statistics or some other
#' statistics used to resume.}}
#'
#' @keywords internal
#' @export
ergm_state <- function(x, ...) UseMethod("ergm_state")

#' @describeIn ergm_state a method for constructing an ergm_state from an [`edgelist`] object and an empty [`network`].
#' @export
ergm_state.edgelist <- function(x, nw0, model=NULL, proposal=NULL, stats=NULL, ...){
  response <- if(ncol(x)>=3) colnames(x)[3]
  out <- list()
  if(is.matrix(x)){
    x <- as_tibble(x)
    colnames(x) <- c(".tail",".head",response)
    x <- as.edgelist(x, n=network.size(nw0), directed=is.directed(nw0), bipartite=nw0%n%"bipartite", loops=has.loops(nw0), vnames=nw0 %v% "vertex.names", output="tibble")
  }
  out$el <- x
  if(!is.null(response)){
    out$el <- out$el[out$el[[response]]!=0,]
    mode(out$el[[3]]) <- "double" # If network is empty, may default to a list().
  }
  out$nw0 <- nw0
  out$nw0[,] <- FALSE
  out$model <- model
  out$proposal <- proposal
  out$stats <- as.double(stats)
  structure(out, class="ergm_state")
}

#' @describeIn ergm_state a method for constructing an ergm_state from a matrix object and an empty [`network`].
#' @export
ergm_state.matrix <- ergm_state.edgelist

#' @describeIn ergm_state a method for constructing an ergm_state from a [`network`] object.
#' @export
ergm_state.network <- function(x, response=NULL, model=NULL, proposal=NULL, stats=NULL, ...){
  NVL(response) <- x %ergmlhs% "response"
  x %ergmlhs% "response" <- response
  el <- as.edgelist(x, attrname=response, output="tibble")
  nw0 <- x
  nw0[,] <- FALSE
  ergm_state(el, nw0, model=model, proposal=proposal, stats=stats, ...)
}

#' @describeIn ergm_state a method for updating an ergm_state.
#' @export
ergm_state.ergm_state <- function(x, el=NULL, nw0=NULL, response=NULL, model=NULL, proposal=NULL, stats=NULL, ...){
  # TODO: Implement sanity checks.
  if(!is.null(nw0)){
    if(is.network(nw0)) x$nw0 <- nw0
    else stop("New nw0 is not a network object.")
  }

  if(!is.null(el)){
    if(is(el, "tibble_edgelist")) x$el <- el
    else stop("New el is not a tibble-style edgelist.")
  }

  if(!is.null(response) && (ncol(x$el)<3 || names(x$el)[3]!=response))
    stop("Attempting to update ergm_state object with a non-matching response attribute.")

  if(!is.null(model)){
    if(is(model, "ergm_model")) x$model <- model
    else stop("New model is not an ergm_model.")
  }

  if(!is.null(proposal)){
    if(is(proposal, "ergm_proposal")) x$proposal <- proposal
    else stop("New proposal is not an ergm_proposal.")
  }

  if(!is.null(stats)) x$stats <- as.double(stats)
  x
}

#' @rdname ergm_state
#' @export
is.ergm_state <- function(x){
  is(x, "ergm_state")
}

#' @rdname ergm_state
#' @export
as.edgelist.ergm_state <- function(x,...){
  x$el
}

#' @rdname ergm_state
#' @export
as.matrix.ergm_state <- function(x,matrix.type=NULL,...){
  matrix.type<-match.arg(matrix.type,c("adjacency","incidence","edgelist"))
  if(matrix.type!="edgelist") stop("ergm_state can only be converted to an edgelist at this time")

  as.matrix(x$el)
}

#' @rdname ergm_state
#' @export
as.network.ergm_state <- function(x, ..., populate=TRUE){
  if(!populate) x$nw0
  else update(x$nw0,x$el)
}

#' @describeIn ergm_state Note that this method fails when
#'   `na.omit=FALSE`, since missing edges are not stored.
#' @param na.omit Whether missing edges should be counted. Note that
#'   missing edge information is not stored.
#' @export
network.edgecount.ergm_state <- function(x, na.omit=TRUE,...){
  if(!na.omit) stop("ergm_state cannot store missing edges.")
  nrow(x$el)
}

#' @describeIn ergm_state Note that this method fails with
#'   its default argument, since missing edges are not stored.
#' @export
network.dyadcount.ergm_state <- function(x, na.omit=TRUE,...){
  if(na.omit) stop("ergm_state cannot store missing edges.")
  network.dyadcount(x$nw0)
}

#' @rdname ergm_state
#' @export
network.size.ergm_state <- function(x,...){
  network.size(x$nw0)
}

#' @describeIn ergm_state A stub that returns 0.
#' @export
network.naedgecount.ergm_state <- function(x,...){
  0
}

#' @rdname ergm_state
#' @export
`%ergmlhs%.ergm_state` <- function(lhs, setting){
  lhs$nw0 %ergmlhs% setting
}

#' @rdname ergm_state
#' @export
`%ergmlhs%<-.ergm_state` <- function(lhs, setting, value){
  lhs$nw0 %ergmlhs% setting <- value
  lhs
}

#' @rdname ergm_state
#' @export
as.rlebdm.ergm_state <- function(x, ...){
  as.rlebdm(x$el, ...)
}

#' @rdname ergm_state
#' @export
as.ergm_model.ergm_state <- function(x, ...) x$model

#' @rdname ergm_state
#' @export
is.curved.ergm_state <- function(object, ...) is.curved(object$model, ...)

#' @rdname ergm_state
#' @export
param_names.ergm_state <- function(object, ...) param_names(object$model, ...)

#' @rdname ergm_state
#' @export
nparam.ergm_state <- function(object, ...) nparam(object$model, ...)
