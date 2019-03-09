#  File R/network.update.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
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

#' @describeIn ergm-deprecated Use [update.network()] instead.
#' @export network.update
network.update<-function(nw, newmatrix, matrix.type=NULL, output="network", ignore.nattr=c("bipartite","directed","hyper","loops","mnext","multiple","n"), ignore.vattr=c()){
  .Deprecate_once("update.network")
  nw[,] <- FALSE

  if(is.null(matrix.type)){
    warning("Don't leave matrix type to chance! Pass matrix.type to network.update!")
    matrix.type <- which.matrix.type(newmatrix)
    if(nrow(newmatrix)==0){matrix.type <- "edgelist"}
  }
  
  if(matrix.type=="adjacency" && all(newmatrix%in%c(0,1))){
    nw[,] <- newmatrix
  }else if(matrix.type=="edgelist" && !is.null(newmatrix) && nrow(newmatrix)>0){
    add.edges(nw,tail=newmatrix[,1],head=newmatrix[,2])
  }
  if(!is.null(output) && output=="edgelist.compressed") 
    {nw <- as.edgelist.compressed(nw)}
  nw
}


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
#'   \code{\link[network]{which.matrix.type}} function.
#' 
#' @param attrname For a network with edge weights gives the name of
#'   the edge attribute whose names to set.
#' 
#' @param ignore.nattr Character vector of the names of network-level
#'   attributes to ignore when updating network objects (defaults to
#'   standard network properties).
#'
#' @param \dots Additional arguments; currently unused.
#' 
#' @param ignore.vattr Character vector of the names of vertex-level
#'   attributes to ignore when updating network objects.
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
update.network <- function(object, new, matrix.type=NULL, attrname=NULL, ..., ignore.nattr=c("bipartite","directed","hyper","loops","mnext","multiple","n"), ignore.vattr=c()){
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
  }else if(matrix.type=="edgelist" && !is.null(new) && nrow(new)>0){
    if(!is.null(attrname)){
      names.eval <- rep(list(attrname), nrow(new))
      vals.eval <- {tmp <- new[,3]; mode(tmp) <- "list"; tmp}
    }else{
      names.eval <- vals.eval <- NULL
    }
    add.edges(object,tail=new[,1],head=new[,2],names.eval=names.eval, vals.eval=vals.eval)
  }
  
  object
}



## FIXME: as.edgelist.compressed and as.network.uncompressed should be
## deleted as soon as network.update() is defunct-ed.
as.edgelist.compressed<-function(x, attrname=NULL, force.bipartite=FALSE, ...){
  .Deprecated(msg="No longer used.")
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


as.network.uncompressed<-function(x, 
        na.rm=FALSE, edge.check=FALSE, ...){
  .Deprecated(msg="No longer used.")
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
