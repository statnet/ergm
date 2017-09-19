#  File R/InitConstraint.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
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

# Baseline constraint incorporating network attributes such as
# directedness, bipartitedness, and self-loops.
InitConstraint..attributes <- function(conlist, lhs.nw, ...){
  conlist$.attributes<-list()

  conlist$.attributes$free.dyads <- function(){
    n <- network.size(lhs.nw)
    ## NB: Free dyad RLE matrix is stored in a column-major order for
    ## consistency with R.
    d <-
      if(has.loops(lhs.nw)) rep(rep(rle(TRUE),n,scale="run"),n,scale="run")
      else do.call(c, lapply(seq_len(n), function(i) rep(rle(c(TRUE,FALSE,TRUE)), c(i-1, 1, n-i),scale="run")))
    
    if(is.bipartite(lhs.nw)){
      n1 <- lhs.nw%n%"bipartite"
      n2 <- n - n1
      
      d <- d &
        c(rep(rep(rle(c(FALSE)),n,scale="run"),n1,scale="run"),
          rep(rep(rle(c(TRUE,FALSE)),c(n1,n2),scale="run"),n2,scale="run"))
    }
    
    if(!is.directed(lhs.nw)){
      d <- d &
        do.call(c, lapply(seq_len(n), function(i) rep(rle(c(TRUE,FALSE)), c(i-1, n-i+1),scale="run")))
    }
    
    rlebdm(compact.rle(d), n)
  }

  conlist
}

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
    as.rlebdm(as.edgelist(is.na(lhs.nw)), matrix.type="edgelist")
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
    n <- network.size(lhs.nw)
    a <- lhs.nw %v% attrname
    a <- rle(a)
    rlebdm(compact.rle(do.call(c,
                        mapply(function(blen,bend){rep(rle(c(FALSE,TRUE,FALSE)), c(bend-blen, blen, n-bend))},
                               a$lengths, cumsum(a$lengths), SIMPLIFY=FALSE)
                        )), n)
  }
  
  conlist
}



InitConstraint.fixedas<-function(conlist, lhs.nw, present=NULL, absent=NULL,...){
	if(is.null(present) & is.null(absent))
		stop(paste("fixedas constraint takes at least one argument, either present or absent or both."), call.=FALSE)
	if(!is.null(present)){
		if(is.network(present)){
			present <- as.edgelist(present)
		}
		if(!is.matrix(present)){
			stop("Argument 'present' in fixedas constraint should be either a network or edgelist")
		}
	}
	if(!is.null(absent)){
		if(is.network(absent)){
			absent <- as.edgelist(absent)
		}
		if(!is.matrix(absent)){
			stop("Argument 'absent' in fixedas constraint should be either a network or edgelist")
		}
	}
	conlist$fixedas$free.dyads<-function(){
          # FixedEdgeList
          fixed <- rbind(present,absent)
          if(any(duplicated(fixed))){
            stop("Dyads cannot be fixed at both present and absent")
          }

          attr(fixed, "n") <- network.size(lhs.nw)
          !as.rlebdm(fixed, matrix.type="edgelist")
        }
        
	conlist
}


InitConstraint.fixallbut<-function(conlist, lhs.nw, free.dyads=NULL,...){
	if(is.null(free.dyads))
		stop(paste("fixallbut constraint takes one required argument free.dyads and one optional argument fixed.state"), call.=FALSE)
	

		if(is.network(free.dyads)){
			free.dyads <- as.edgelist(free.dyads)
		}
		
		if(!is.matrix(free.dyads)){
			stop("Argument 'free.dyads' in fixallbut constraint should be either a network or edgelist")
		}
	
#	
#	if(!is.null(fixed.state)){
#		if(length(fixed.state)==1)
#			rep(fixed.state,nrow(fixed.dyads))
#		if(length(fixed.state != nrow(fixed.dayds)))
#			stop("fixed.state should be a vector of length equals to the nubmer of dyads in fixed.dyads")
#		if(!all(fixed.state %in% c(0,1)))
#			stop("fixed.state should be a vector of 0,1")
#	}
#	
#	
	conlist$fixallbut$free.dyads<-function(){ 
          attr(free.dyads, "n") <- network.size(lhs.nw)
          as.rlebdm(free.dyads, matrix.type="edgelist")
	}
	conlist
}





#ergm.ConstraintImplications("edges", c())
