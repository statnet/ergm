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
#' # Draw a random directed network with 25 nodes
#' g <- network(25)
#'
#' # Draw a random undirected network with density 0.1
#' g <- network(25, directed=FALSE, density=0.1)
#'
#' # Draw a random bipartite network with 4 actors and 6 events and density 0.1
#' g <- network(10, bipartite=4, directed=FALSE, density=0.1)
#'
#' # Draw a random directed network with 25 nodes and 50 edges
#' g <- network(25, numedges=50)
#' @importFrom network as.network
#' @export
as.network.numeric<-function(x,
    directed = TRUE,
    hyper = FALSE, loops = FALSE, multiple = FALSE, bipartite = FALSE,
    ignore.eval = TRUE, names.eval = NULL,
    edge.check = FALSE,
    density=NULL, init=NULL, numedges=NULL, ...){
  # Producing an informative error for each of the following invalid or unsupported inputs
  if(loops || multiple || hyper)
    stop("Generating multigraphs, or networks with self-loops or hyperedges is not supported at this time.")
  if(NVL3(density, .<0 || .>1, FALSE))
    stop("Density of graph cannot be either negative or greater then 1")
  if(NVL3(numedges, round(.)!=., FALSE))
    stop("The number of edges must be an integer")
  ## # TODO: After network() with match.call() is on CRAN, enable this.
  ## if(!missing(directed) && directed && bipartite != FALSE)
  ##   stop("Generating directed bipartite networks is not supported at this time.")
  ## if(bipartite!=FALSE && missing(directed)){
  if(bipartite!=FALSE && directed){
    directed <- FALSE
    warning_once("Bipartite network specified: assuming undirected. Pass ", sQuote("directed=FALSE")," to silence this warning. This behavior may change in the future.")
  }
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

  if(ndyads > 2^(53-4))
    stop("The number of possible edges cannot be greater than 2^49.")

  if(NVL(numedges, 0) > ndyads)
    stop("The number of edges cannot be greater than the number of possible edges.")

  if(missing(density)){
    if(missing(init)){
      #     So the expected number of ties is the same as
      #     the number of nodes
      density <- nb1/ndyads
    }else{
      density <- exp(init)/(1+exp(init))
    }
  }

  if(is.null(numedges))
    numedges <- rbinom(1,ndyads,density)

  index <- sample.int(ndyads, numedges)

  if(directed){
    tails <- (index %/% (nb2-1L)) + 1L - (index%%(nb2-1)==0L)
    heads <- index - (tails-1L)*(nb2-1L)
    heads <- heads + (heads>=tails)
  }else if(bipartite){
    heads <- index %% nb2
    heads[heads==0L] <- nb2
    tails <- 1L + ((index - heads)/nb2)
    heads <- heads + nb1
  }else{
    difvi <- ceiling(sqrt(8L*(ndyads - index)+9L))
    tails <- nb2 - (difvi + (difvi%%2L==0L) - 1L)/2L
    heads <- index + tails*(tails+1L)/2L - (tails-1L)*nb2
  }

  el <- structure(cbind(as.integer(tails), as.integer(heads)),
                  n = as.integer(x), bipartite = if(bipartite) as.integer(bipartite))
  as.network(el, directed = directed, matrix.type = "edgelist")
}
