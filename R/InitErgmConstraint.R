#  File R/InitErgmConstraint.R in package ergm, part of the Statnet suite
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
# constraint lists (each prependend with "InitErgmConstraint")
#         <edges>                   <odegreedist>
#         <degrees>=<nodedegrees>   <bd>
#         <degreesTetrad>           <idegrees>
#         <degreesHexad>            <odegrees>
#         <degreedist>              <hamming>
#         <idegreedist>            <observed>
#============================================================================

##########################################################################################
# Each of the <InitErgmConstraint.X> functions accepts an existing constraint list, 'conlist',
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

#' Get multiple vertex attributes at once and paste them together.
#'
#' @param x,na.omit,null.na,unlist see [network::get.vertex.attribute()].
#' @param attrnames a character vector of one or more vertex attribute
#'   names; see [network::get.vertex.attribute()].
#' @param sep an optional character vector of length 1 to use as a
#'   separator for attribute values.
#'
#' @return If `sep` is `NULL`, a list with an element for each element
#'   of `attrnames` containing the vertex attribute vector for that
#'   attribute. Otherwise, if a character vector, convert the
#'   attribute values to strings and join them with `sep` as the
#'   separator.
get.vertex.attributes <- function(x, attrnames, na.omit = FALSE, null.na = TRUE, unlist=TRUE, sep=NULL){
  a <- lapply(attrnames, get.vertex.attribute, x=x, na.omit=na.omit, null.na=null.na, unlist=unlist)
  if(!is.null(sep)) a <- do.call(paste, c(a, sep=sep))
  a
}

# Baseline constraint incorporating network attributes such as
# directedness, bipartitedness, and self-loops.
InitErgmConstraint..attributes <- function(lhs.nw, ...){
  list(
    free_dyads = {
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
    },
    constrain = character(0),
    dependence = FALSE)
}

InitErgmConstraint.edges<-function(lhs.nw, ...){
   if(length(list(...)))
     stop(paste("Edge count constraint does not take arguments at this time."), call.=FALSE)
   list(dependence = TRUE)
}
#ergm.ConstraintImplications("edges", c())

InitErgmConstraint.degrees<-InitErgmConstraint.nodedegrees<-function(lhs.nw, ...){
   if(length(list(...)))
     stop(paste("Vertex degrees constraint does not take arguments at this time."), call.=FALSE)
   list(dependence = TRUE, constrain = "degrees")
}

#ergm.ConstraintImplications("degrees", c("edges", "idegrees", "odegrees", "idegreedist", "odegreedist", "degreedist", "bd"))

InitErgmConstraint.odegrees<-function(lhs.nw, ...){
   if(length(list(...)))
     stop(paste("Vertex odegrees constraint does not take arguments at this time."), call.=FALSE)
   if(!is.directed(lhs.nw)) stop("Vertex odegrees constraint is only meaningful for directed networks.", call.=FALSE)
   list(dependence = TRUE)
}
#ergm.ConstraintImplications("odegrees", c("edges", "odegreedist"))

InitErgmConstraint.idegrees<-function(lhs.nw, ...){
   if(length(list(...)))
     stop(paste("Vertex idegrees constraint does not take arguments at this time."), call.=FALSE)
   if(!is.directed(lhs.nw)) stop("Vertex idegrees constraint is only meaningful for directed networks.", call.=FALSE)
   list(dependence = TRUE)
}
#ergm.ConstraintImplications("idegrees", c("edges", "idegreedist"))

InitErgmConstraint.b1degrees<-function(lhs.nw, ...){
   if(length(list(...)))
     stop(paste("B1 vertex degrees constraint does not take arguments at this time."), call.=FALSE)
   if(!is.bipartite(lhs.nw)) stop("B1 vertex degrees constraint is only meaningful for bipartite networks.", call.=FALSE)
   list(dependence = TRUE)
}
#ergm.ConstraintImplications("b1degrees", c("edges"))

InitErgmConstraint.b2degrees<-function(lhs.nw, ...){
   if(length(list(...)))
     stop(paste("B2 vertex degrees constraint does not take arguments at this time."), call.=FALSE)
   if(!is.bipartite(lhs.nw)) stop("B2 vertex degrees constraint is only meaningful for bipartite networks.", call.=FALSE)
   list(dependence = TRUE)
}
#ergm.ConstraintImplications("b2degrees", c("edges"))

InitErgmConstraint.degreedist<-function(lhs.nw, ...){
   if(length(list(...)))
     stop(paste("Degree distribution constraint does not take arguments at this time."), call.=FALSE)
   list(dependence = TRUE)
}
#ergm.ConstraintImplications("degreedist", c("edges", "idegreedist", "odegreedist"))


InitErgmConstraint.idegreedist<-function(lhs.nw, ...){
   if(length(list(...)))
     stop(paste("InDegree distribution constraint does not take arguments at this time."), call.=FALSE)
   list(dependence = TRUE)
}
#ergm.ConstraintImplications("idegreedist", c("edges"))


InitErgmConstraint.odegreedist<-function(lhs.nw, ...){
   if(length(list(...)))
     stop(paste("OutDegree distribution constraint does not take arguments at this time."), call.=FALSE)
   list(dependence = TRUE)
}
#ergm.ConstraintImplications("odegreedist", c("edges"))


InitErgmConstraint.bd<-function(lhs.nw, attribs=NULL, maxout=NA, maxin=NA, minout=NA, minin=NA, ...){
   if(nargs()>6)
     stop(paste("Bounded degrees constraint takes at most 5 arguments; ",nargs()-1," given.",sep=""), call.=FALSE)
   list(attribs=attribs,maxout=maxout,maxin=maxin,minout=minout,minin=minin)
}
#ergm.ConstraintImplications("bd", c())

InitErgmConstraint.hamming<-function(lhs.nw, ...){
   if(length(list(...)))
     stop(paste("Hamming distance constraint does not take arguments at this time."), call.=FALSE)
   list(dependence = TRUE)
}
#ergm.ConstraintImplications("hamming", c())

InitErgmConstraint.observed <- function(lhs.nw, ...){
  if(length(list(...)))
    stop(paste("Toggle non-observed constraint does not take arguments at this time."), call.=FALSE)
  list(free_dyads = as.rlebdm(as.edgelist(is.na(lhs.nw))),
       dependence = FALSE)
}
#ergm.ConstraintImplications("observed", c())

InitErgmConstraint.blockdiag<-function(lhs.nw, attrname=NULL, ...){
  if(length(list(...)))
    stop(paste("Block diagonal constraint takes one argument at this time."), call.=FALSE)
  list(attrname=attrname,
       free_dyads = {
         n <- network.size(lhs.nw)
         a <- lhs.nw %v% attrname
         if(NVL(lhs.nw%n%"bipartite",0)){
           bip <- lhs.nw %n% "bipartite"
           ea <- a[seq_len(bip)]
           aa <- a[bip+seq_len(n-bip)]
           if(length(rle(ea)$lengths)!=length(unique(rle(ea)$values)) || length(rle(aa)$lengths)!=length(unique(rle(aa)$values))) stop("Current implementation of block-diagonal sampling requires that the blocks of the egos and the alters be contiguous. See help('ergm-constraints') for more information.")
           
           tmp <- .double.rle(ea, aa)
           el <- tmp$lengths1
           al <- tmp$lengths2
           
           o <- rlebdm(c(rep(rle(FALSE), bip*n, scale="run"),
                         do.call(c,rep(
                                     mapply(function(blen,bend){rep(rle(c(FALSE,TRUE,FALSE)), c(bend-blen, blen, n-bend), scale="run")},
                                            el, cumsum(el), SIMPLIFY=FALSE),
                                     al)
                                 )), n)
           # Future-proofing: in case it's bipartite directed, add
           # both thte blocks and their transposes. (If undirected,
           # it'll get filtered out by the .attributes constraints.)
           ot <- rlebdm(c(do.call(c,rep(
                                      mapply(function(blen,bend){rep(rle(c(FALSE,TRUE,FALSE)), c(bip+bend-blen, blen, n-bip-bend), scale="run")},
                                             al, cumsum(al), SIMPLIFY=FALSE),
                                      el)
                                  ),
                          rep(rle(FALSE), (n-bip)*n, scale="run")), n)
           o | ot
         }else{
           a <- rle(a)
           rlebdm(compact.rle(do.call(c,rep(
                                          mapply(function(blen,bend){rep(rle(c(FALSE,TRUE,FALSE)), c(bend-blen, blen, n-bend), scale="run")},
                                                 a$lengths, cumsum(a$lengths), SIMPLIFY=FALSE),
                                          a$lengths)
                                      )), n)
         }
       },
       dependence = FALSE)
}



InitErgmConstraint.fixedas<-function(lhs.nw, present=NULL, absent=NULL,...){
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
  list(
    free_dyads = {
      # FixedEdgeList
      fixed <- as.edgelist(rbind(present,absent),
                           n=lhs.nw%n%"n",
                           directed=lhs.nw%n%"directed",
                           bipartite=lhs.nw%n%"bipartite",
                           loops=lhs.nw%n%"loops")
      if(any(duplicated(fixed))){
        stop("Dyads cannot be fixed at both present and absent")
      }

      !as.rlebdm(fixed)
    },
    dependence = FALSE)
}


InitErgmConstraint.fixallbut<-function(lhs.nw, free.dyads=NULL,...){
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
  list(
    free_dyads = { 
      fixed <- as.edgelist(free.dyads,
                           n=lhs.nw%n%"n",
                           directed=lhs.nw%n%"directed",
                           bipartite=lhs.nw%n%"bipartite",
                           loops=lhs.nw%n%"loops")
      as.rlebdm(free.dyads)
    },
    dependence = FALSE)
}


InitErgmConstraint.dyadnoise<-function(lhs.nw, p01, p10, ...){
  if(length(list(...)))
    stop(paste("Dyadic noise \"constraint\" takes one argument at this time."), call.=FALSE)

  if(((length(p01) != 1 || length(p10) != 1) &&
        any(dim(as.matrix(lhs.nw, matrix.type="adjacency")) != c(dim(p01),dim(p10)))))
    stop("p01 and p10 must be either scalars or matrices of the same dimension as the adjacency matrices of the LHS network.")
  
  list(p01=p01, p10=p10)
}



#ergm.ConstraintImplications("edges", c())


InitErgmConstraint.egocentric <- function(lhs.nw, attrname=NULL, direction = c("both", "out", "in")){
  direction <- match.arg(direction)
  if(!is.directed(lhs.nw) && direction!="both"){
    stop("Directed egocentric constraint cannot be used for an undirected network.")
  }
  n <- network.size(lhs.nw)
  a <- ( # Are that node's dyads toggleable?
    if(is.null(attrname)) get.vertex.attribute(lhs.nw, "na")
    else !get.vertex.attribute(lhs.nw, attrname)
  )

  list(
    free_dyads = {
      # Remember: column-major order.

      rlea <- rle(a)
      
      fd <- rlebdm(switch(direction,
                          `out` = rep(rlea, n),
                          `in` = rep(rlea, rep(n, length(rlea$lengths)), scale="run"),
                          `both` = rep(rlea, n) & rep(rlea, rep(n, length(rlea$lengths)), scale="run")),
                   n)
    },
    dependence = FALSE
  )
}
