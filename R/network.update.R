#  File R/network.update.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
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

network.update<-function(nw, newmatrix, matrix.type=NULL, output="network")
{
#  print(paste("old:",network.edgecount(nw)," new:", nrow(newmatrix),collapse=" "))
  unw <- network.copy(as.network.uncompressed(nw))
  if(is.null(matrix.type)){
    warning("Don't leave matrix type to chance! Pass matrix.type to network.update!")
    matrix.type <- which.matrix.type(newmatrix)
    if(nrow(newmatrix)==0){matrix.type <- "edgelist"}
  }
  if(matrix.type=="adjacency" && all(newmatrix%in%c(0,1))){
    unw[,] <- newmatrix
  }else if(matrix.type=="edgelist"){
#  cnw <- as.matrix.network(nw,matrix.type="edgelist")
#  unw[cnw[,2],cnw[,1]] <- 0
#  unw[,] <- 0
#  eid<-vector()
#  for(i in 1:network.size(nw)){  
#    eid <- c(eid,get.edgeIDs(unw,i))
#  }
   eid <- c(unlist(unw$iel),unlist(unw$oel))
   delete.edges(unw,eid)
   if(!is.null(newmatrix) && nrow(newmatrix)>0){
#   unw[newmatrix] <- 1
     
    # *** don't forget - edges are given as tails -> heads now
    add.edges(unw,tail=newmatrix[,1],head=newmatrix[,2])
   }
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
