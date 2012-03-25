#============================================================================
# This file contains the following 12 functions for initializing empty
# constraint lists (each prependend with "InitConstraint")
#         <edges>                   <outdegreedist>
#         <degrees>=<nodedegrees>   <bd>
#         <degreesTetrad>           <indegrees>
#         <degreesHexad>            <outdegrees>
#         <degreedist>              <hamming>
#         <indegreedist>            <observed>
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
#ergm.ConstraintImplications("degrees", c("edges", "indegrees", "outdegrees", "indegreedist", "outdegreedist", "degreedist", "bd"))

InitConstraint.outdegrees<-InitConstraint.outdegrees<-function(conlist, lhs.nw, ...){
   if(length(list(...)))
     stop(paste("Vertex outdegrees constraint does not take arguments at this time."), call.=FALSE)
   if(!is.directed(lhs.nw)) stop("Vertex outdegrees constraint is only meaningful for directed networks.", call.=FALSE)
   conlist$outdegrees<-list()
   conlist
}
#ergm.ConstraintImplications("outdegrees", c("edges", "outdegreedist"))

InitConstraint.indegrees<-InitConstraint.indegrees<-function(conlist, lhs.nw, ...){
   if(length(list(...)))
     stop(paste("Vertex indegrees constraint does not take arguments at this time."), call.=FALSE)
   if(!is.directed(lhs.nw)) stop("Vertex indegrees constraint is only meaningful for directed networks.", call.=FALSE)
   conlist$indegrees<-list()
   conlist
}
#ergm.ConstraintImplications("indegrees", c("edges", "indegreedist"))

InitConstraint.b1degrees<-InitConstraint.b1degrees<-function(conlist, lhs.nw, ...){
   if(length(list(...)))
     stop(paste("B1 vertex degrees constraint does not take arguments at this time."), call.=FALSE)
   if(!is.bipartite(lhs.nw)) stop("B1 vertex degrees constraint is only meaningful for bipartite networks.", call.=FALSE)
   conlist$b1degrees<-list()
   conlist
}
#ergm.ConstraintImplications("b1degrees", c("edges"))

InitConstraint.b2degrees<-InitConstraint.b2degrees<-function(conlist, lhs.nw, ...){
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
#ergm.ConstraintImplications("degreedist", c("edges", "indegreedist", "outdegreedist"))


InitConstraint.indegreedist<-function(conlist, lhs.nw, ...){
   if(length(list(...)))
     stop(paste("InDegree distribution constraint does not take arguments at this time."), call.=FALSE)
   conlist$indegreedist<-list()
   conlist
}
#ergm.ConstraintImplications("indegreedist", c("edges"))


InitConstraint.outdegreedist<-function(conlist, lhs.nw, ...){
   if(length(list(...)))
     stop(paste("OutDegree distribution constraint does not take arguments at this time."), call.=FALSE)
   conlist$outdegreedist<-list()
   conlist
}
#ergm.ConstraintImplications("outdegreedist", c("edges"))


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
