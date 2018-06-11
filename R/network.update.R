#  File R/network.update.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2018 Statnet Commons
#######################################################################
#============================================================================
# This file contains the following 3 functions used to update networks:
#          <network.update>
#          <as.edgelist.compressed>
#          <as.network.uncompressed>
#===========================================================================



###############################################################################
# The <network.update> function returns the given network with with only the
# ties specified by a given matrix
#
# --PARAMETERS--
#   nw         : a network object
#   newmatrix  : the matrix specifying the new set of ties with which to
#                update 'nw' 
#   matrix.type: the type of matrix that 'newmatrix' is, as "adjacency" or
#                "edgelist"; default=which.matrix.type(newmatrix)
#   output     : a string indicating whether the output should be an
#                edgelist (using "edgelist.compressed") or should be a 
#                network (using any other string); default="network"
#
# --RETURNED--
#   unw:  the updated network, having only those ties specified by 'newmatrix'
#
###############################################################################

#' Create an empty copy of a network object
#' 
#' Initializes an empty network with the same vertex and network
#' attributes as the original network, but no edges.
#'
#' @param nw a [`network`] object
#' @param ignore.nattr character vector of the names of network-level
#'   attributes to ignore when updating network objects (defaults to
#'   standard network properties)
#' @param ignore.vattr character vector of the names of vertex-level
#'   attributes to ignore when updating network objects
empty_network <- function(nw, ignore.nattr=c("bipartite","directed","hyper","loops","mnext","multiple","n"), ignore.vattr=c()){
  if(network.edgecount(nw)==0) return(nw)
  
  unw <- network.initialize(n=network.size(nw), directed = is.directed(nw), hyper = is.hyper(nw), loops = has.loops(nw),
         multiple = is.multiplex(nw), bipartite = nw %n% "bipartite")
  for(a in setdiff(list.network.attributes(nw),ignore.nattr)) unw <- set.network.attribute(unw, a, get.network.attribute(nw, a, unlist=FALSE))
  for(a in setdiff(list.vertex.attributes(nw),ignore.vattr)) unw <- set.vertex.attribute(unw, a, get.vertex.attribute(nw, a, unlist=FALSE))
  unw
}

#' Replace the sociomatrix in a network object
#' 
#' Replaces the edges in a network object with the edges corresponding
#' to the sociomatrix specified by \code{newmatrix}.  See
#' \code{\link{ergm}} for more information.
#' 
#' 
#' @param nw a \code{\link[network]{network}} object. See
#'   documentation for the \code{\link[network]{network}} package.
#' @param newmatrix Either an adjacency matrix (a matrix of zeros and
#'   ones indicating the presence of a tie from i to j) or an edgelist
#'   (a two-column matrix listing origin and destination node numbers
#'   for each edge; note that in an undirected matrix, the first
#'   column should be the smaller of the two numbers).
#' @param matrix.type One of "adjacency" or "edgelist" telling which
#'   type of matrix \code{newmatrix} is.  Default is to use the
#'   \code{\link[network]{which.matrix.type}} function.
#' @param output Currently unused.
#' @param ignore.nattr character vector of the names of network-level
#'   attributes to ignore when updating network objects (defaults to
#'   standard network properties)
#' @param ignore.vattr character vector of the names of vertex-level
#'   attributes to ignore when updating network objects
#' @return \code{\link{network.update}} returns a new
#'   \code{\link[network]{network}} object with the edges specified by
#'   \code{newmatrix} and network and vertex attributes copied from
#'   the input network \code{nw}. Input network is not modified.
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
#' network.update(flomarriage, rand.mat, matrix.type="adjacency")
#' # Try this with an edgelist
#' rand.mat <- as.matrix.network.edgelist(flomarriage)[1:5,]
#' network.update(flomarriage, rand.mat, matrix.type="edgelist")
#' 
#' @export network.update
network.update<-function(nw, newmatrix, matrix.type=NULL, output="network", ignore.nattr=c("bipartite","directed","hyper","loops","mnext","multiple","n"), ignore.vattr=c()){
  unw <- empty_network(nw, ignore.nattr=ignore.nattr, ignore.vattr=ignore.vattr)

  if(is.null(matrix.type)){
    warning("Don't leave matrix type to chance! Pass matrix.type to network.update!")
    matrix.type <- which.matrix.type(newmatrix)
    if(nrow(newmatrix)==0){matrix.type <- "edgelist"}
  }
  
  if(matrix.type=="adjacency" && all(newmatrix%in%c(0,1))){
    unw[,] <- newmatrix
  }else if(matrix.type=="edgelist" && !is.null(newmatrix) && nrow(newmatrix)>0){
    add.edges(unw,tail=newmatrix[,1],head=newmatrix[,2])
  }
  if(!is.null(output) && output=="edgelist.compressed") 
    {unw <- as.edgelist.compressed(unw)}
  unw
}


###############################################################################
# The <as.edgelist.compressed> function converts a network 'x' into the edgelist
# 'out' described below; this is a copy of <as.edgelist.san>
#
# --PARAMETERS--
#   x              : a network object, or a list of such
#   attrname       : optionally, the name of an edge attribute to use for edge
#                    values; default=NULL
#   force.bipartite: whether ?? if 'x' is not already bipartite(T or F); default=FALSE; if TRUE,
#                    this appears to merely create the 'input must be a network'
#                    warning, before finishing up as if this were FALSE
#
# --RETURNED--
#   out: x, as an edgelist with attributes for
#       n                : the network size
#       directed         : whether the network is directed (T or F)
#       vnames           : the vertex names
#       vertex.attributes: a list of the vertex attributes
#       bipartite        : whether the network is bipartite (T or F)
#
###############################################################################

as.edgelist.compressed<-function(x, attrname=NULL, force.bipartite=FALSE, ...){
  #In case of lists, process independently
  if(is.list(x) && !inherits(x,"network"))
    return(lapply(x,as.edgelist.compressed, attrname=attrname, force.bipartite=force.bipartite))
  #Begin with network objects
  if(inherits(x,"network")){
    out<-as.matrix.network.edgelist(x,attrname=attrname)
#   if(!is.directed(x)){
#    out <- out[1:(nrow(x)/2),]
#   }
    if(NCOL(out)==2)                        #If needed, add edge values
      out<-cbind(out,rep(1,NROW(out)))
    attr(out,"n")<-network.size(x)
    attr(out,"directed")<-is.directed(x)
    attr(out,"vnames")<-network.vertex.names(x)
    van<-list.vertex.attributes(x)
    if(length(van)>0){
     va <- vector(mode = "list", length(van))
     for (i in (1:length(van))){ 
      va[[i]]<-get.vertex.attribute(x,van[i],unlist=TRUE)
     }
     names(va)<-van
     attr(out,"vertex.attributes")<-va
    }
    if(is.bipartite(x))
      attr(out,"bipartite")<-get.network.attribute(x,"bipartite")
    else if(force.bipartite)
      out<-as.edgelist.compressed(out,attrname=attrname,force.bipartite=force.bipartite)
  }else{
    warning("as.edgelist.compressed input must be network, or list thereof.\n Returning the original object.\n")
    return(x)
  }
  #Return the result
  out
}



###############################################################################
# The <as.network.uncompressed> function is basically the inverse of the above
# <as.edgelist.compressed> function
#
# --PARAMETERS--
#   x         : a compressed network or a network
#   edge.check: whether computationally expensive checks of the legality
#               of submitted edges should be performed (T or F); default=FALSE
#
# --IGNORED PARAMTERS--
#   na.rm:  whether NA valuse should be removed for ??; default=FALSE
#   ...  :  additional parameters for flexibility
#
# --RETURNED--
#   x: the original network if it is already uncompressed or if 'x' is neither
#      a compressed or uncompressed network
#   g: the uncompressed version of x
#
###############################################################################

as.network.uncompressed<-function(x, 
        na.rm=FALSE, edge.check=FALSE, ...){
  #Initialize the network object
  if(inherits(x,"network")){return(x)}
  if(is.null(attr(x,"vnames"))){
   warning("as.network.uncompressed input must be a compressed network, or a network.\n Returning the original object.\n")
   return(x)
  }
  n<-attr(x,"n")
  directed<-attr(x,"directed")
  g<-network.initialize(n,directed=directed)
  #Call the specific coercion routine, depending on matrix type
# g<-network.edgelist(x,g,na.rm=na.rm,edge.check=edge.check)
  g<-add.edges(g,as.list(x[,1]),as.list(x[,2]),edge.check=edge.check)
  va <- attr(x,"vertex.attributes")
  if(length(va)>0){
   for (i in (1:length(va))){ 
    g <- set.vertex.attribute(g,names(va)[i], va[[i]])
   }
  }
  #Return the result
  g
}
