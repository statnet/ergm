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


InitConstraint.bd<-function(conlist, lhs.nw, attribs=0, maxout=0, maxin=0, minout=0, minin=0){
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
  conlist
}
#ergm.ConstraintImplications("observed", c())

InitConstraint.ranks<-function(conlist, lhs.nw, ...){
   if(length(list(...)))
     stop(paste("Rank constraint does not take arguments at this time."), call.=FALSE)
   conlist$ranks<-list()
   conlist
}
