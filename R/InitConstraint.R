#  File ergm/R/InitConstraint.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
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

InitConstraint.degreessimple<-function(conlist){
   if (nargs()>1)
     stop(paste("Vertex degrees constraint does not take arguments at this time."), call.=FALSE)
   conlist$degreessimple<-list()
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

InitConstraint.atleast<-function(conlist,nw=NULL){
  if(is.null(nw)) stop("Formation constraint ``atleast'' requires a baseline network.",call.=FALSE)
  if(network.naedgecount(nw)) stop("Baseline network passed to formation constraint ``atleast'' may not have missing dyads.")
  conlist$atleast<-list(nw=nw)
  conlist
}

InitConstraint.atmost<-function(conlist,nw=NULL){
  if(is.null(nw)) stop("Dissolution constraint ``atmost'' requires a baseline network.",call.=FALSE)
  if(network.naedgecount(nw)) stop("Baseline network passed to dissolution constraint ``atleast'' may not have missing dyads.")
  conlist$atmost<-list(nw=nw)
  conlist
}
