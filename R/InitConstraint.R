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


## List of which constraints make which constraints redundant.
ConstraintImplications<-list(edges=character(0),
                             degrees=c("edges","indegrees","outdegrees","degreedist","bd"),
                             degreedist=c("edges","indegrees"),
                             indegreedist=c("edges"),
                             outdegreedist=c("edges"),
                             bd=character(0),
                             indegrees=c("edges"),
                             outdegrees=c("edges"),
                             hamming=character(0),
                             observed=character(0))



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

InitConstraint.edges<-function(conlist){
   if (nargs()>1)
     stop(paste("Edge count constraint does not take arguments at this time."), call.=FALSE)
   conlist$edges<-list()
   conlist
}


InitConstraint.degrees<-InitConstraint.nodedegrees<-function(conlist){
   if (nargs()>1)
     stop(paste("Vertex degrees constraint does not take arguments at this time."), call.=FALSE)
   conlist$degrees<-list()
   conlist
}



InitConstraint.degreesTetrad<-function(conlist){
   if (nargs()>1)
     stop(paste("Vertex degreesTetrad constraint does not take arguments at this time."), call.=FALSE)
   conlist$degreesTetrad<-list()
   conlist
}



InitConstraint.degreesHexad<-function(conlist){
   if (nargs()>1)
     stop(paste("Vertex degreesHexad constraint does not take arguments at this time."), call.=FALSE)
   conlist$Hexad<-list()
   conlist
}



InitConstraint.degreedist<-function(conlist){
   if (nargs()>1)
     stop(paste("Degree distribution constraint does not take arguments at this time."), call.=FALSE)
   conlist$degreedist<-list()
   conlist
}



InitConstraint.indegreedist<-function(conlist){
   if (nargs()>1)
     stop(paste("InDegree distribution constraint does not take arguments at this time."), call.=FALSE)
   conlist$indegreedist<-list()
   conlist
}



InitConstraint.outdegreedist<-function(conlist){
   if (nargs()>1)
     stop(paste("OutDegree distribution constraint does not take arguments at this time."), call.=FALSE)
   conlist$outdegreedist<-list()
   conlist
}



InitConstraint.bd<-function(conlist, attribs=0, maxout=0, maxin=0, minout=0, minin=0){
   if (nargs()>6)
     stop(paste("Bounded degrees constraint takes at most 5 arguments; ",nargs()-1," given.",sep=""), call.=FALSE)
   conlist$bd<-list(attribs=attribs,maxout=maxout,maxin=maxin,minout=minout,minin=minin)
   conlist
}



#InitConstraint.indegrees<-function(conlist){
#   if (nargs()>1)
#     stop(paste("Vertex indegrees constraint does not take arguments at this time."), call.=FALSE)
#   conlist$indegrees<-list()
#   conlist
#}



#InitConstraint.outdegrees<-function(conlist){
#   if (nargs()>1)
#     stop(paste("Vertex outdegrees constraint does not take arguments at this time."), call.=FALSE)
#   conlist$outdegrees<-list()
#   conlist
#}



InitConstraint.hamming<-function(conlist){
   if (nargs()>1)
     stop(paste("Hamming distance constraint does not take arguments at this time."), call.=FALSE)
   conlist$hamming<-list()
   conlist
}



InitConstraint.observed <- function(conlist){
  if(nargs()>1)
    stop(paste("Toggle non-observed constraint does not take arguments at this time."), call.=FALSE)
  conlist$observed<-list()
  conlist
}
