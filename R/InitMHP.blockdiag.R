#  File R/InitMHP.blockdiag.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
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



#### WARNING: The following functions also have a copy in tergm. Fixes
#### should be applied to both (for now.)
## FIXME: There is almost certainly a better way to do this.
## TODO: Document functions and export them, for use by tergm.
.consensus.order <- function(x1, x2){
  o <- intersect(x1, x2)
  if(!all(x1[x1 %in% o] == x2[x2 %in% o])) stop("Current implementation of block-diagonal sampling requires the common blocks of egos and blocks of alters to have the same order. See help('ergm-constraints') for more information.")
  o1 <- c(0, which(x1 %in% o),length(x1)+1)
  o2 <- c(0, which(x2 %in% o),length(x2)+1)
  n <- length(o1) - 1
  v <- c()

  sr <- function(from,to){from + seq_len(to-from + 1) - 1}
    
  for(i in seq_len(n)){
    v <- c(v, x1[sr(o1[i]+1,o1[i+1]-1)])
    v <- c(v, x2[sr(o2[i]+1,o2[i+1]-1)])
    v <- na.omit(c(v, x1[o1[i+1]]))
  }
  as.vector(v)
}

.double.rle <- function(a1, a2){
  e1 <- rle(a1)
  e2 <- rle(a2)

  o <- .consensus.order(e1$values, e2$values)

  l1 <- e1$lengths[match(o, e1$values)]
  l1[is.na(l1)] <- 0
  l2 <- e2$lengths[match(o, e2$values)]
  l2[is.na(l2)] <- 0
  
  list(values=o, lengths1=l1, lengths2=l2)
}

#' @title Compute and serialize information needed by the block-diagonal
#' Metropolis-Hastings samplers.
#'
#' @description Given a nodal attribute vector specifying the blocks
#'   returns a list containing information needed to efficiently
#'   sample dyads in those blocks.
#'
#' @param nw the network of interest.
#' @param a either a vector (of any mode) or a vertex attribute of
#'   `nw` whose \eqn{i}th element identifies the block to which
#'   vertex \eqn{i} belongs. Blocks must be continguous.
#' 
#' @return A list with the following elements:
#' \describe{
#'
#' \item{`nblk`}{number of blocks}
#' 
#' \item{`pos` (for unipartite), `b1pos` and `b2pos` (for
#' bipartite)}{vectors of length `nblk+1` such that
#' `(pos[i]+1):(pos[i+1])` are the numbers of the vertices in block
#' `i`.}
#'
#' \item{`cumwt`}{vector of length `nblk` with values such that
#' generating \eqn{U\sim \mathrm{Uniform}(0,1)}{U ~ Uniform(0,1)}
#' selecting the index of the first element in `cumwt` that is greater
#' than \eqn{U} will sample from the blocks in proportion to the
#' number of dyads in the block.}
#'
#' }
#'
#' In addition, an attribute `"ndyads"` is attached, containing the
#' total number of dyads in all blocks put together.
#'
#' `unlist(pack.BlockDiagSampInfo_as_num(nw, blkattr)` can be
#' passed as input to an MH proposal function that knows how to parse
#' it.
#' 
#' @export
pack.BlockDiagSampInfo_as_num <- function(nw, a){
  bip <- nw %n% "bipartite"
  a <- if(length(a)==network.size(nw)) a else get.vertex.attributes(nw, a, sep="\t")
  if(bip){
    ea <- a[seq_len(bip)]
    aa <- a[bip+seq_len(network.size(nw)-bip)]
  
    ## rle() returns contigous runs of values.
    # If we have more runs than unique values, the blocks must not be all contiguous.
    if(length(rle(ea)$lengths)!=length(unique(rle(ea)$values)) || length(rle(aa)$lengths)!=length(unique(rle(aa)$values))) stop("Current implementation of block-diagonal sampling requires that the blocks of the egos and the alters be contiguous. See help('ergm-constraionts') for more information.")

    tmp <- .double.rle(ea, aa)

    eb <- cumsum(c(0,tmp$lengths1)) # upper bounds of ego blocks
    ab <- cumsum(c(0,tmp$lengths2))+bip # upper bounds of alter blocks
    w <- cumsum(tmp$lengths1*tmp$lengths2) # cumulative block weights ~ # dyads in the block
    w <- w/max(w)
    # Note that this automagically takes care of singleton blocks by giving them weight 0.
    out <- list(nblk=length(w),  b1pos=eb, b2pos=ab, cumwt=w)
    attr(out, "ndyads") <- sum(tmp$lengths1*tmp$lengths2)
  }else{
    ## rle() returns contigous runs of values.
    a <- rle(a)
    # If we have more runs than unique values, the blocks must not be all contiguous.
    if(length(a$lengths)!=length(unique(a$values))) stop("Current implementation of block-diagonal sampling requires that the blocks be contiguous.")
    b <- cumsum(c(0,a$lengths)) # upper bounds of blocks
    w <- cumsum(a$lengths*(a$lengths-1)) # cumulative block weights ~ # dyads in the block
    w <- w/max(w)
    # Note that this automagically takes care of singleton blocks by giving them weight 0.
    out <- list(nblk=length(w),  pos=b, cumwt=w)
    attr(out, "ndyads") <- sum(a$lengths*(a$lengths-1)/(if(is.directed(nw)) 1 else 2))
  }
  return(out)
}

InitMHP.blockdiag <- function(arguments, nw){
  BDI <- pack.BlockDiagSampInfo_as_num(nw, arguments$constraints$blockdiag$attrname)  
  list(name = "blockdiag", inputs=unlist(BDI))
}

InitMHP.blockdiagTNT <- function(arguments, nw){
  el <- as.edgelist(nw)
  a <- get.vertex.attributes(nw, arguments$constraints$blockdiag$attrname, sep="\t")
  
  if(any(a[el[,1]]!=a[el[,2]])) stop("Block-diagonal TNT sampler implementation does not support sampling networks with off-block-diagonal ties at this time.")

  BDI <- pack.BlockDiagSampInfo_as_num(nw, arguments$constraints$blockdiag$attrname)
  
  list(name = "blockdiagTNT", inputs=c(attr(BDI,"ndyads"), unlist(BDI)))
}

## Helper function, since the following two have the same body except for the MH_ function.
.InitMHP.blockdiagNonObserved <- function(arguments, nw, ...){
  ## Bipartite is handled seamlessly.
  
  a <- get.vertex.attributes(nw, arguments$constraints$blockdiag$attrname, sep="\t")

  el <- as.edgelist(is.na(nw))
  el <- el[a[el[,1]]==a[el[,2]],,drop=FALSE]
  
  c(list(inputs=ergm.Cprepare.el(el)), list(...))
}


InitMHP.blockdiagNonObserved <- function(arguments, nw){
  .InitMHP.blockdiagNonObserved(arguments, nw, name = "randomtoggleList")
}

InitMHP.blockdiagNonObservedTNT <- function(arguments, nw){
  .InitMHP.blockdiagNonObserved(arguments, nw, name = "listTNT")
}
