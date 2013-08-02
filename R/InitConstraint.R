#  File R/InitConstraint.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
#============================================================================
# This file contains the following 12 functions for initializing empty
# constraint lists (each prependend with "InitConstraint")
#         <edges>                   <odegreedist>
#         <degrees>=<nodedegrees>   <bd>
#         <degreesTetrad>           <idegrees>
#         <degreesHexad>            <odegrees>
#         <degreedist>              <hamming>
#         <idegreedist>            <observed>
#============================================================================

##########################################################################################
# Each of the <InitConstraint.X> functions accepts an existing constraint list, 'conlist',
# and to this adds an empty constraint list for term X; if any arguments are passed besides
# 'conlist", execution will halt.
#
# --PARAMETERS--
#   conlist: a list, presumably of constraints for other terms
#
# --RETURNED--
#   conlist: updated to include the initialized empty constraint list for term X
#
##########################################################################################

InitConstraint.edges<-function(conlist, lhs.nw, ...){
   if(length(list(...)))
     stop(paste("Edge count constraint does not take arguments at this time."), call.=FALSE)
   conlist$edges<-list()
   conlist
}
#ergm.ConstraintImplications("edges", c())

InitConstraint.degrees<-InitConstraint.nodedegrees<-function(conlist, lhs.nw, ...){
   if(length(list(...)))
     stop(paste("Vertex degrees constraint does not take arguments at this time."), call.=FALSE)
   conlist$degrees<-list()
   conlist
}

#ergm.ConstraintImplications("degrees", c("edges", "idegrees", "odegrees", "idegreedist", "odegreedist", "degreedist", "bd"))

InitConstraint.odegrees<-function(conlist, lhs.nw, ...){
   if(length(list(...)))
     stop(paste("Vertex odegrees constraint does not take arguments at this time."), call.=FALSE)
   if(!is.directed(lhs.nw)) stop("Vertex odegrees constraint is only meaningful for directed networks.", call.=FALSE)
   conlist$odegrees<-list()
   conlist
}
#ergm.ConstraintImplications("odegrees", c("edges", "odegreedist"))

InitConstraint.idegrees<-function(conlist, lhs.nw, ...){
   if(length(list(...)))
     stop(paste("Vertex idegrees constraint does not take arguments at this time."), call.=FALSE)
   if(!is.directed(lhs.nw)) stop("Vertex idegrees constraint is only meaningful for directed networks.", call.=FALSE)
   conlist$idegrees<-list()
   conlist
}
#ergm.ConstraintImplications("idegrees", c("edges", "idegreedist"))

InitConstraint.b1degrees<-function(conlist, lhs.nw, ...){
   if(length(list(...)))
     stop(paste("B1 vertex degrees constraint does not take arguments at this time."), call.=FALSE)
   if(!is.bipartite(lhs.nw)) stop("B1 vertex degrees constraint is only meaningful for bipartite networks.", call.=FALSE)
   conlist$b1degrees<-list()
   conlist
}
#ergm.ConstraintImplications("b1degrees", c("edges"))

InitConstraint.b2degrees<-function(conlist, lhs.nw, ...){
   if(length(list(...)))
     stop(paste("B2 vertex degrees constraint does not take arguments at this time."), call.=FALSE)
   if(!is.bipartite(lhs.nw)) stop("B2 vertex degrees constraint is only meaningful for bipartite networks.", call.=FALSE)
   conlist$b2degrees<-list()
   conlist
}
#ergm.ConstraintImplications("b2degrees", c("edges"))

InitConstraint.degreedist<-function(conlist, lhs.nw, ...){
   if(length(list(...)))
     stop(paste("Degree distribution constraint does not take arguments at this time."), call.=FALSE)
   conlist$degreedist<-list()
   conlist
}
#ergm.ConstraintImplications("degreedist", c("edges", "idegreedist", "odegreedist"))


InitConstraint.idegreedist<-function(conlist, lhs.nw, ...){
   if(length(list(...)))
     stop(paste("InDegree distribution constraint does not take arguments at this time."), call.=FALSE)
   conlist$idegreedist<-list()
   conlist
}
#ergm.ConstraintImplications("idegreedist", c("edges"))


InitConstraint.odegreedist<-function(conlist, lhs.nw, ...){
   if(length(list(...)))
     stop(paste("OutDegree distribution constraint does not take arguments at this time."), call.=FALSE)
   conlist$odegreedist<-list()
   conlist
}
#ergm.ConstraintImplications("odegreedist", c("edges"))


InitConstraint.bd<-function(conlist, lhs.nw, attribs=NULL, maxout=NA, maxin=NA, minout=NA, minin=NA){
   if(nargs()>6)
     stop(paste("Bounded degrees constraint takes at most 5 arguments; ",nargs()-1," given.",sep=""), call.=FALSE)
   conlist$bd<-list(attribs=attribs,maxout=maxout,maxin=maxin,minout=minout,minin=minin)
   conlist
}
#ergm.ConstraintImplications("bd", c())

InitConstraint.hamming<-function(conlist, lhs.nw, ...){
   if(length(list(...)))
     stop(paste("Hamming distance constraint does not take arguments at this time."), call.=FALSE)
   conlist$hamming<-list()
   conlist
}
#ergm.ConstraintImplications("hamming", c())


InitConstraint.observed <- function(conlist, lhs.nw, ...){
  if(length(list(...)))
    stop(paste("Toggle non-observed constraint does not take arguments at this time."), call.=FALSE)
  conlist$observed<-list()

  conlist$observed$free.dyads <- function(){
    standardize.network(is.na(lhs.nw))
  }
  conlist
}
#ergm.ConstraintImplications("observed", c())

InitConstraint.blockdiag<-function(conlist, lhs.nw, attrname=NULL, ...){
  if(length(list(...)))
    stop(paste("Block diagonal constraint takes one argument at this time."), call.=FALSE)
  conlist$blockdiag <- list(attrname=attrname)
  
  # This definition should "remember" attrname and lhs.nw.
  conlist$blockdiag$free.dyads <- function(){
    a <- lhs.nw %v% attrname
    el <- do.call(rbind,tapply(seq_along(a),INDEX=list(a),simplify=FALSE,FUN=function(i) do.call(rbind,lapply(i,function(j) cbind(j,i)))))
    el <- el[el[,1]!=el[,2],]
    el <- as.edgelist(el, n=network.size(lhs.nw), directed=is.directed(lhs.nw))
    # standardize.network() not needed here, since el is already in standard order.
    network.update(lhs.nw, el, matrix.type="edgelist")
  }
  
  conlist
}
