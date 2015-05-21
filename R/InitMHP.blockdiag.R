#  File R/InitMHP.blockdiag.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2015 Statnet Commons
#######################################################################
########################################################################
# Each of the <InitMHP.X> functions initializes and returns a
# MHproposal list; when appropriate, proposal types are checked against
# covariates and network types for 1 of 2 side effects: to print warning
# messages or to halt execution (only <InitMHP.nobetweengroupties> can
# halt execution)
#
# --PARAMETERS--
#   arguments: is ignored by all but <InitMHP.nobetweengroupties>,
#              where 'arguments' is used to get the nodal attributes
#              via <get.node.attr>
#   nw       : the network given by the model
#
# --RETURNED--
#   MHproposal: a list containing:
#        name   : the name of the proposal
#        inputs : a vector to be passed to the proposal
#        package: is "ergm"
#
############################################################################

InitMHP.blockdiag <- function(arguments, nw){
  if(is.bipartite(nw)) stop("Block-diagonal sampling is not implemented for bipartite networks at this time.")
  # rle() returns contigous runs of values.
  a <- rle(nw %v% arguments$constraints$blockdiag$attrname)
  # If we have more runs than unique values, the blocks must not be all contiguous.
  if(length(a$lengths)!=length(unique(a$values))) stop("Current implementation of block-diagonal sampling requires that the blocks be contiguous.")
  b <- cumsum(c(0,a$lengths)) # upper bounds of blocks
  w <- cumsum(a$lengths*(a$lengths-1)) # cumulative block weights ~ # dyads in the block
  w <- w/max(w)
  # Note that this automagically takes care of singleton blocks by giving them weight 0.
  
  MHproposal <- list(name = "blockdiag", inputs=c(length(b)-1,  b, w))
  MHproposal
}

InitMHP.blockdiagTNT <- function(arguments, nw){
  if(is.bipartite(nw)) stop("Block-diagonal sampling is not implemented for bipartite networks at this time.")

  el <- as.edgelist(nw)
  a <- nw %v% arguments$constraints$blockdiag$attrname
  
  if(any(a[el[,1]]!=a[el[,2]])) stop("Block-diagonal TNT sampler implementation does not support sampling networks with off-block-diagonal ties at this time.")

  
  # rle() returns contigous runs of values.
  a <- rle(a)
  # If we have more runs than unique values, the blocks must not be all contiguous.
  if(length(a$lengths)!=length(unique(a$values))) stop("Current implementation of block-diagonal sampling requires that the blocks be contiguous.")

  nd <- sum(a$lengths*(a$lengths-1)/(if(is.directed(nw)) 1 else 2))
  b <- cumsum(c(0,a$lengths)) # upper bounds of blocks
  w <- cumsum(a$lengths*(a$lengths-1)) # cumulative block weights ~ # dyads in the block
  w <- w/max(w)
  # Note that this automagically takes care of singleton blocks by giving them weight 0.
  
  MHproposal <- list(name = "blockdiagTNT", inputs=c(nd, length(b)-1,  b, w))
  MHproposal
}

InitMHP.blockdiagNonObserved <- function(arguments, nw){
  if(is.bipartite(nw)) stop("Block-diagonal sampling is not implemented for bipartite networks at this time.")
  # rle() returns contigous runs of values.
  a <- nw %v% arguments$constraints$blockdiag$attrname
  # If we have more runs than unique values, the blocks must not be all contiguous.
  if(length(rle(a)$lengths)!=length(unique(rle(a)$values))) stop("Current implementation of block-diagonal sampling requires that the blocks be contiguous.")

  el <- as.edgelist(is.na(nw))

  el <- el[a[el[,1]]==a[el[,2]],,drop=FALSE]
  
  MHproposal <- list(name = "randomtoggleList", inputs=ergm.Cprepare.el(el))
  MHproposal
}

InitMHP.blockdiagNonObservedTNT <- function(arguments, nw){
  if(is.bipartite(nw)) stop("Block-diagonal sampling is not implemented for bipartite networks at this time.")
  # rle() returns contigous runs of values.
  a <- nw %v% arguments$constraints$blockdiag$attrname
  # If we have more runs than unique values, the blocks must not be all contiguous.
  if(length(rle(a)$lengths)!=length(unique(rle(a)$values))) stop("Current implementation of block-diagonal sampling requires that the blocks be contiguous.")

  el <- as.edgelist(is.na(nw))

  el <- el[a[el[,1]]==a[el[,2]],,drop=FALSE]
  
  MHproposal <- list(name = "listTNT", inputs=ergm.Cprepare.el(el))
  MHproposal
}
