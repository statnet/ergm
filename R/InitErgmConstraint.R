#  File R/InitErgmConstraint.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
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
    free_dyads = function(){
      n <- network.size(lhs.nw)
      storage.mode(n) <- "integer"
      ## NB: Free dyad RLE matrix is stored in a column-major order for
      ## consistency with R.
      d <-
        if(is.directed(lhs.nw)){
          if(has.loops(lhs.nw)){
            compress(structure(list(lengths=rep(n,n), values=rep(TRUE,n)), class="rle"))
          }else{
            structure(list(lengths=c(1L,rep(c(n,1L),n-1L)), values=c(rep(c(FALSE, TRUE),n-1L),FALSE)), class="rle")
          }
        }else if(is.bipartite(lhs.nw)){
          b1 <- as.integer(lhs.nw%n%"bipartite")
          b2 <- n - b1
          compress(structure(list(lengths=c(rep(n,b1), rep(c(b1,b2),b2)), values=c(rep(FALSE, b1), rep(c(TRUE,FALSE),b2))),class="rle"))
        }else{
          if(has.loops(lhs.nw)){
            vals <- c(rep(c(TRUE,FALSE),n-1L),TRUE)
            lens <- integer(2L*(n-1L)+1L)
            for(i in seq_len(n-1L)){
              lens[2L*i-1L] <- i
              lens[2L*i] <- n-i
            }
            lens[2L*n-1L] <- n
          }else{
            vals <- c(rep(c(FALSE,TRUE),n-1L),FALSE)
            lens <- integer(2L*(n-1L)+1L)
            for(i in seq_len(n-1L)){
              lens[2L*i-1L] <- n-i+1L
              lens[2L*i] <- i
            }
            lens[2L*n-1L] <- 1L
          }          
          structure(list(lengths=lens,values=vals), class="rle")
        }
      rlebdm(d, n)
    },
    implies = ".attributes",
    dependence = FALSE)
}

InitErgmConstraint.edges<-function(lhs.nw, ...){
   if(...length())
     ergm_Init_abort(paste("Edge count constraint does not take arguments at this time."))
   list(dependence = TRUE, implies = "edges")
}

InitErgmConstraint.degrees<-InitErgmConstraint.nodedegrees<-function(lhs.nw, ...){
   if(...length())
     ergm_Init_abort(paste("Vertex degrees constraint does not take arguments at this time."))
   list(dependence = TRUE, constrain = "degrees", implies = c("degrees", "edges", "idegrees", "odegrees", "idegreedist", "odegreedist", "degreedist", "bd"))
}

InitErgmConstraint.odegrees<-function(lhs.nw, ...){
   if(...length())
     ergm_Init_abort(paste("Vertex odegrees constraint does not take arguments at this time."))
   if(!is.directed(lhs.nw)) ergm_Init_abort("Vertex odegrees constraint is only meaningful for directed networks.")
   list(dependence = TRUE, implies = c("odegrees", "edges", "odegreedist"))
}

InitErgmConstraint.idegrees<-function(lhs.nw, ...){
   if(...length())
     ergm_Init_abort(paste("Vertex idegrees constraint does not take arguments at this time."))
   if(!is.directed(lhs.nw)) ergm_Init_abort("Vertex idegrees constraint is only meaningful for directed networks.")
   list(dependence = TRUE, implies = c("idegrees", "edges", "idegreedist"))
}

InitErgmConstraint.b1degrees<-function(lhs.nw, ...){
   if(...length())
     ergm_Init_abort(paste("B1 vertex degrees constraint does not take arguments at this time."))
   if(!is.bipartite(lhs.nw)) ergm_Init_abort("B1 vertex degrees constraint is only meaningful for bipartite networks.")
   list(dependence = TRUE, implies = c("b1degrees", "edges"))
}

InitErgmConstraint.b2degrees<-function(lhs.nw, ...){
   if(...length())
     ergm_Init_abort(paste("B2 vertex degrees constraint does not take arguments at this time."))
   if(!is.bipartite(lhs.nw)) ergm_Init_abort("B2 vertex degrees constraint is only meaningful for bipartite networks.")
   list(dependence = TRUE, implies = c("b2degrees", "edges"))
}

InitErgmConstraint.degreedist<-function(lhs.nw, ...){
   if(...length())
     ergm_Init_abort(paste("Degree distribution constraint does not take arguments at this time."))
   list(dependence = TRUE, implies = c("degreedist", "edges", "idegreedist", "odegreedist"))
}

InitErgmConstraint.idegreedist<-function(lhs.nw, ...){
   if(...length())
     ergm_Init_abort(paste("InDegree distribution constraint does not take arguments at this time."))
   list(dependence = TRUE, implies = c("idegreedist", "edges"))
}

InitErgmConstraint.odegreedist<-function(lhs.nw, ...){
   if(...length())
     ergm_Init_abort(paste("OutDegree distribution constraint does not take arguments at this time."))
   list(dependence = TRUE, implies = c("odegreedist", "edges"))
}

InitErgmConstraint.bd<-function(lhs.nw, attribs=NULL, maxout=NA, maxin=NA, minout=NA, minin=NA, ...){
   if(nargs()>6)
     ergm_Init_abort(paste("Bounded degrees constraint takes at most 5 arguments; ",nargs()-1," given.",sep=""))
   if(...length()) ergm_Init_abort(paste0("Unrecognised argument(s) ", paste.and(names(list(...)), oq="'", cq="'"),"."))

   if(!is.directed(lhs.nw) && is.null(attribs) && length(maxout) == 1 && is.na(maxin) && is.na(minout) && is.na(minin)) {
     constrain <- "bdmax"
   } else {
     constrain <- "bd"
   }

   list(constrain=constrain, attribs=attribs, maxout=maxout, maxin=maxin, minout=minout, minin=minin)
}

InitErgmConstraint.blocks <- function(lhs.nw, attr = NULL, fmat = NULL, ...) {
  if(...length()) ergm_Init_abort(paste0("Unrecognised argument(s) ", paste.and(names(list(...)), oq = "'", cq = "'"), ".")) 

  if(is.directed(lhs.nw)) {
    ergm_Init_abort("blocks constraint does not currently support directed networks.")
  }

  ## dyad-independent constraint
  dependence <- FALSE    
  
  free_dyads <- function() {
    n <- as.integer(network.size(lhs.nw))
    
    if(is.bipartite(lhs.nw)) {
      b1 <- as.integer(lhs.nw %n% "bipartite")
      b2 <- n - b1
      
      bd_row_nodecov <- NVL2(attr, ergm_get_vattr(attr, lhs.nw, bip = "b1"), integer(b1))
      bd_col_nodecov <- NVL2(attr, ergm_get_vattr(attr, lhs.nw, bip = "b2"), integer(b2))    
    } else {
      b1 <- 0L
      b2 <- 0L
      
      bd_row_nodecov <- NVL2(attr, ergm_get_vattr(attr, lhs.nw), integer(n))
      bd_col_nodecov <- bd_row_nodecov    
    }
    
    bd_row_levels <- sort(unique(bd_row_nodecov))
    bd_col_levels <- sort(unique(bd_col_nodecov))
      
    bd_row_nodecov <- match(bd_row_nodecov, bd_row_levels)
    bd_col_nodecov <- match(bd_col_nodecov, bd_col_levels)
    
    fmat <- NVL(fmat, matrix(FALSE, nrow = length(bd_row_levels), ncol = length(bd_col_levels)))
    amat <- matrix(!fmat, nrow = nrow(fmat), ncol = ncol(fmat))
    
    rle_list <- vector(mode = "list", length = length(bd_col_levels))
    for(i in seq_along(bd_col_levels)) {
      rle_list[[i]] <- rle(c(amat[bd_row_nodecov,i], logical(b2)))
    }
    
    lens <- vector(mode = "list", length = n)
    vals <- vector(mode = "list", length = n)

    for(i in seq_len(n)) {
      if(i > b1) {
        lens[[i]] <- rle_list[[bd_col_nodecov[i - b1]]]$lengths
        vals[[i]] <- rle_list[[bd_col_nodecov[i - b1]]]$values
      } else {
        lens[[i]] <- n
        vals[[i]] <- FALSE    
      }
    }

    rlebdm(compress(structure(list(lengths=unlist(lens), values=unlist(vals)), class="rle")), n)
  }
  
  list(dependence = dependence, free_dyads = free_dyads, attr = attr, fmat = fmat)  
}


InitErgmConstraint.hamming<-function(lhs.nw, ...){
   if(...length())
     ergm_Init_abort(paste("Hamming distance constraint does not take arguments at this time."))
   list(dependence = TRUE)
}

InitErgmConstraint.observed <- function(lhs.nw, ...){
  if(...length())
    ergm_Init_abort(paste("Toggle non-observed constraint does not take arguments at this time."))
  list(free_dyads = as.rlebdm(as.edgelist(is.na(lhs.nw))),
       dependence = FALSE, implies = c("observed"))
}

InitErgmConstraint.fixedas<-function(lhs.nw, present=NULL, absent=NULL,...){
  if(is.null(present) && is.null(absent))
    ergm_Init_abort(paste("fixedas constraint takes at least one argument, either present or absent or both."))
  if(!is.null(present) && !is.network(present) && !is.matrix(present))
    ergm_Init_abort("Argument 'present' in fixedas constraint should be either a network or edgelist")
  if(!is.null(absent) && !is.network(present) && !is.matrix(absent))
    ergm_Init_abort("Argument 'absent' in fixedas constraint should be either a network or edgelist")

  list(
    free_dyads = function(){
      if(is.network(present)) present <- as.edgelist(present)
      if(is.network(absent)) absent <- as.edgelist(absent)

      # FixedEdgeList
      fixed <- as.edgelist(rbind(present,absent),
                           n=lhs.nw%n%"n",
                           directed=lhs.nw%n%"directed",
                           bipartite=lhs.nw%n%"bipartite",
                           loops=lhs.nw%n%"loops")
      if(any(duplicated(fixed))){
        ergm_Init_abort("Dyads cannot be fixed at both present and absent")
      }

      !as.rlebdm(fixed)
    },
    dependence = FALSE)
}


InitErgmConstraint.fixallbut<-function(lhs.nw, free.dyads=NULL,...){
  if(is.null(free.dyads))
    ergm_Init_abort(paste("fixallbut constraint takes one required argument free.dyads"))
  if(!is.network(free.dyads) && !is.matrix(free.dyads))
    ergm_Init_abort("Argument 'free.dyads' in fixallbut constraint should be either a network or edgelist")
  
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
    free_dyads = function(){
      if(is.network(free.dyads)) free.dyads <- as.edgelist(free.dyads)
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
  if(...length())
    stop(paste("Dyadic noise \"constraint\" takes one argument at this time."), call.=FALSE)

  if(((length(p01) != 1 || length(p10) != 1) &&
        any(dim(as.matrix(lhs.nw, matrix.type="adjacency")) != c(dim(p01),dim(p10)))))
    stop("p01 and p10 must be either scalars or matrices of the same dimension as the adjacency matrices of the LHS network.")
  
  list(p01=p01, p10=p10)
}

InitErgmConstraint.egocentric <- function(lhs.nw, attr=NULL, direction = c("both", "out", "in")){
  direction <- match.arg(direction)
  if(!is.directed(lhs.nw) && direction!="both")
    stop("Directed egocentric constraint cannot be used for an undirected network.")

  list(
    free_dyads = function(){
      n <- network.size(lhs.nw)
      a <- ( # Are that node's dyads toggleable?
        if(is.null(attr)) get.vertex.attribute(lhs.nw, "na")
        else !as.vector(ergm_get_vattr(attr, lhs.nw, accept="logical"))
      )

      # Remember: column-major order.

      rlea <- rle(a)
      
      fd <- rlebdm(switch(direction,
                          `out` = rep(rlea, n),
                          `in` = rep(rlea, rep(n, length(rlea$lengths)), scale="run"),
                          `both` = compress(rep(rlea, n) & rep(rlea, rep(n, length(rlea$lengths)), scale="run"))), # The others are already compressed by rep().
                   n)
    },
    dependence = FALSE
  )
}

InitErgmConstraint.Dyads<-function(lhs.nw, fix=NULL, vary=NULL,...){
  if(is.null(fix) & is.null(vary))
    ergm_Init_abort(paste("Dyads constraint takes at least one argument, either",sQuote("fix"),"or",sQuote("vary"),"or both."))

  list(
    free_dyads = function(){
      fd <- lapply(list(fix=fix,vary=vary),
                   function(f){
                     if(!is.null(f)){
                       f[[3]] <- f[[2]]
                       f[[2]] <- lhs.nw
                       m <- ergmMPLE(f, expand.bipartite=TRUE, output="array")$predictor
                       m <- m!=0
                       m[is.na(m)] <- FALSE
                       if(!is.directed(lhs.nw)){
                         m <- m | aperm(m, c(2L,1L,3L))
                       }
                       lapply(seq_len(dim(m)[3]), function(i) as.rlebdm(m[,,i]))
                     }
                   })
      fd$fix <- if(length(fd$fix)) fd$fix %>% map(`!`) %>% reduce(`&`)
      fd$vary <- if(length(fd$vary)) fd$vary %>% reduce(`|`)
      fd <- Reduce(`|`, fd)

      compress(fd)
    },
    dependence = FALSE
  )
}
