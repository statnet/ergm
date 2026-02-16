#  File R/network.update.R in package ergm, part of the Statnet suite of
#  packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################

#' Update the edges in a network based on a matrix
#' 
#' Replaces the edges in a [`network`] object with the edges corresponding
#' to the sociomatrix or edge list specified by \code{new}.
#' 
#' 
#' @param object a [`network`] object.
#' 
#' @param new Either an adjacency matrix (a matrix of values
#'   indicating the presence and/or the value of a tie from i to j) or
#'   an edge list (a two-column matrix listing origin and destination
#'   node numbers for each edge, with an optional third column for the
#'   value of the edge).
#' 
#' @param matrix.type One of `"adjacency"` or `"edgelist"` telling
#'   which type of matrix \code{new} is.  Default is to use the
#'   [which.matrix.type()] function.
#' 
#' @param attrname For a network with edge weights gives the name of
#'   the edge attribute whose names to set.
#'
#' @param \dots Additional arguments; currently unused.
#' 
#' @return A new [`network`] object with the edges specified by
#'   \code{new} and network and vertex attributes copied from
#'   the input network `object`. Input network is not modified.
#' 
#' @seealso [ergm()], [`network`]
#' @keywords models
#' @examples
#' 
#' #
#' data(florentine)
#' #
#' # test the network.update function
#' #
#' # Create a Bernoulli network
#' rand.net <- network(network.size(flomarriage))
#' # store the sociomatrix 
#' rand.mat <- rand.net[,]
#' # Update the network
#' update(flomarriage, rand.mat, matrix.type="adjacency")
#' # Try this with an edgelist
#' rand.mat <- as.matrix.network.edgelist(flomarriage)[1:5,]
#' update(flomarriage, rand.mat, matrix.type="edgelist")
#' 
#' @export
update.network <- function(object, ...){
  update_network(object, ...)
}

#' @describeIn update.network dispatcher for network update based on the type of updating information.
#' @export
update_network <- function(object, new, ...){
  UseMethod("update_network", new)
}

#' @describeIn update.network a method for updating a network based on a matrix-form edgelist
#' @export
update_network.matrix_edgelist <- function(object, new, attrname=if(ncol(new)>2) names(new)[3], ...){
  NVL(colnames(new)) <- c(".tail",".head",attrname)
  new <- as_tibble(new)
  update_network(object, new, attrname=attrname, ...)
}

#' @describeIn update.network a method for updating a network based on an edgelist
#' @export
update_network.data.frame <- function(object, new, attrname=if(ncol(new)>2) names(new)[3], ...){
  # Empty the network.
  object[,] <- FALSE
  if(nrow(new)){
    if(!is.null(attrname)){
      names.eval <- rep(list(attrname), nrow(new))
      vals.eval <- {tmp <- new[[3]]; mode(tmp) <- "list"; tmp}
    }else{
      names.eval <- vals.eval <- NULL
    }
    add.edges(object,tail=new[[1]],head=new[[2]],names.eval=names.eval, vals.eval=vals.eval)
  }
  object
}

#' @describeIn update.network a method for updating a network based on a matrix
#' @export
update_network.matrix <- function(object, new, matrix.type=NULL, attrname=NULL, ...){
  if(is.null(matrix.type)){
    warning("Don't leave matrix type to chance! Pass matrix.type to update.network!")
    matrix.type <- which.matrix.type(new)
    if(nrow(new)==0){matrix.type <- "edgelist"}
  }

  if(! matrix.type%in%c("adjacency","edgelist")) stop("Only edge lists and adjacency matrices are supporeted at this time.")

  # Empty the network.
  object[,] <- FALSE

  if(matrix.type=="adjacency"){
    object[,,names.eval=attrname,add.edges=TRUE] <- new
  }else{
    object <- update_network.matrix_edgelist(object, new, matrix.type=matrix.type, attrname=attrname, ...)
  }
  
  object
}

#' @describeIn update.network a method for updating a network based on an [`ergm_state`] object.
#' @export
update_network.ergm_state <- function(object, new, ...){
  update.network(object, new$el)  
  object
}
