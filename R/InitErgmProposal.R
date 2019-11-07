#  File R/InitErgmProposal.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################
#===========================================================================
# The <InitErgmProposal> file contains the following 24 functions for
# initializing the proposal object; each is prepended with 'InitErgmProposal.'
#       <randomtoggle>      <CondOutDegreeDist> 
#       <TNT>               <ConstantEdges>    
#       <CondInDegree>      <CondDegree>
#       <CondOutDegree>     <HammingTNT>   
#       <CondDegreeTetrad>         <HammingConstantEdges>
#       <CondDegreeHexad>            <randomtoggleNonObserved>
#       <CondDegreeDist>          <nobetweengroupties>
#       <CondInDegreeDist>  
#============================================================================


########################################################################
# Each of the <InitErgmProposal.X> functions initializes and returns a
# proposal list; when appropriate, proposal types are checked against
# covariates and network types for 1 of 2 side effects: to print warning
# messages or to halt execution (only <InitErgmProposal.nobetweengroupties> can
# halt execution)
#
# --PARAMETERS--
#   arguments: is ignored by all but <InitErgmProposal.nobetweengroupties>,
#              where 'arguments' is used to get the nodal attributes
#              via <get.node.attr>
#   nw       : the network given by the model
#
# --RETURNED--
#   proposal: a list containing:
#        name   : the name of the proposal
#        inputs : a vector to be passed to the proposal
#        package: is "ergm"
#
############################################################################
InitErgmProposal.randomtoggle <- function(arguments, nw) {
  proposal <- list(name = "randomtoggle", inputs=NULL)
  proposal
}

InitErgmProposal.TNT <- function(arguments, nw) {
  proposal <- list(name = "TNT", inputs=NULL)
  proposal
}

InitErgmProposal.StratTNT <- function(arguments, nw) {
  if(is.null(arguments$attr))
    ergm_Init_abort("The ", sQuote("attr"), " argument to ", sQuote("StratTNT"), " is required (and must be named).")

  pmat <- arguments$pmat

  if(!is.bipartite(nw)) {
    nodecov <- ergm_get_vattr(arguments$attr, nw)
    levels <- sort(unique(nodecov))
    nodecov <- match(nodecov, levels)

    # default is matrix of 1s
    if(is.null(pmat)) pmat <- matrix(1, nrow = length(levels), ncol = length(levels))
    
    if(!is.matrix(pmat) || !is.numeric(pmat)) ergm_Init_abort("The ", sQuote("pmat"), " argument to ", sQuote("StratTNT"), " must be a numeric matrix.")
    
    if(NROW(pmat) != length(levels) || NCOL(pmat) != length(levels))
      ergm_Init_abort("For unipartite networks, the ", sQuote("pmat"), " argument to ", sQuote("StratTNT"), " must be a square matrix with number of rows and number of columns both equal to the number of unique values of the ", sQuote("attr"), " argument.")
  
    # if undirected unipartite, then symmetrize and set the sub-diagonal to zero
    if(!is.directed(nw)) {
      pmat <- (pmat + t(pmat))/2
      pmat[!upper.tri(pmat, diag=TRUE)] <- 0
    }
  } else {
    b1nodecov <- ergm_get_vattr(arguments$attr, nw, bip="b1")
    b2nodecov <- ergm_get_vattr(arguments$attr, nw, bip="b2")
    
    b1levels <- sort(unique(b1nodecov))
    b2levels <- sort(unique(b2nodecov))
    
    # default is matrix of 1s    
    if(is.null(pmat)) pmat <- matrix(1, nrow = length(b1levels), ncol = length(b2levels))    
    
    if(!is.matrix(pmat) || !is.numeric(pmat)) ergm_Init_abort("The ", sQuote("pmat"), " argument to ", sQuote("StratTNT"), " must be a numeric matrix.")    
    
    if(NROW(pmat) != length(b1levels) || NCOL(pmat) != length(b2levels))
      ergm_Init_abort("For bipartite networks, the ", sQuote("pmat"), " argument to ", sQuote("StratTNT"), " must be a matrix with number of rows equal to the number of unique values of the ", sQuote("attr"), " argument on the first bipartition, and number of columns equal to the number of unique values of the ", sQuote("attr"), " argument on the second bipartition.")
      
    # shift b2 codes so there is no overlap with b1 codes
    nodecov <- c(match(b1nodecov, b1levels), match(b2nodecov, b2levels) + length(b1levels))
  }

  if(any(is.na(pmat)) || any(pmat < 0) || all(pmat <= 0)) {
    ergm_Init_abort("The ", sQuote("pmat"), " argument to ", sQuote("StratTNT"), " must be a numeric matrix with non-negative entries, at least one of which is strictly positive.  It cannot contain any missing data.")
  }
  
  # renormalize to probability matrix
  pmat <- pmat/sum(pmat)

  # record the tail and head attr code for each mixing type with positive probability
  tailattrs <- NULL
  headattrs <- NULL
  probvec <- NULL
  
  for(i in 1:NROW(pmat)) {
    for(j in 1:NCOL(pmat)) {
      if(pmat[i,j] > 0) {
        tailattrs <- c(tailattrs, i)
        headattrs <- c(headattrs, j)
        probvec <- c(probvec, pmat[i,j])
      }
    }
  }

  # if bipartite, shift b2 codes (which always correspond to heads) so there is no overlap with b1 codes
  if(is.bipartite(nw)) {
    headattrs <- headattrs + length(b1levels)
  }

  # order mixing types so largest probabilities occur first (hopefully save a little time in the C code)
  o <- rev(order(probvec))
  
  tailattrs <- tailattrs[o]
  headattrs <- headattrs[o]
  probvec <- cumsum(probvec[o])

  # record the number of mixing types and the number of unique attr codes
  nmixingtypes <- length(probvec)
  ncodes <- max(nodecov)
  
  # record the attr codes in the order of the nodal indices
  codesbynodeindex <- nodecov
  
  # record the nodal indices grouped (and counted) by attr code
  nodeindicesbycode <- NULL
  nodecountsbycode <- NULL
  
  for(i in 1:ncodes) {
    w <- which(nodecov == i)
    nodeindicesbycode <- c(nodeindicesbycode, w)
    nodecountsbycode <- c(nodecountsbycode, length(w))
  }
  
  # check that all mixing types with positive probability have at least 1 dyad; error if not
  for(i in 1:nmixingtypes) {
    if(nodecountsbycode[tailattrs[i]] == 0 || nodecountsbycode[headattrs[i]] == 0 || (tailattrs[i] == headattrs[i] && nodecountsbycode[tailattrs[i]] == 1)) {
      ergm_Init_abort("Mixing types with positive proposal probability must have at least one dyad.")
    }
  }
  
  # awkwardly force everything into one big vector for the C code...
  inputs <- c(nmixingtypes, tailattrs, headattrs, probvec, ncodes, nodecountsbycode, nodeindicesbycode, codesbynodeindex)
    
  proposal <- list(name = "StratTNT", inputs=inputs)
  proposal
}

InitErgmProposal.CondDegree <- function(arguments, nw) {
  proposal <- list(name = "CondDegree", inputs=NULL)
  proposal
}
InitErgmProposal.CondDegreeMix <- function(arguments, nw) {
  proposal <- list(name = "CondDegreeMix",
    inputs=get.vertex.attribute(nw,arguments$constraints$degreesmix$attrib))
  proposal
}

InitErgmProposal.CondOutDegree <- function(arguments, nw) {
  proposal <- list(name = "CondOutDegree", inputs=NULL)
  if (!is.directed(nw)) # Really, this should never trigger, since the InitErgmConstraint function should check.
    ergm_Init_abort("The CondOutDegree proposal function does not work with an",
          "undirected network.")
  
  proposal
}

InitErgmProposal.CondInDegree <- function(arguments, nw) {
  proposal <- list(name = "CondInDegree", inputs=NULL)
  if (!is.directed(nw)) # Really, this should never trigger, since the InitErgmConstraint function should check.
    ergm_Init_abort("The CondInDegree proposal function does not work with an",
          "undirected network.")
  proposal
}

InitErgmProposal.CondB1Degree <- function(arguments, nw) {
  proposal <- list(name = "CondB1Degree", inputs=NULL)
  if (!is.bipartite(nw)) # Really, this should never trigger, since the InitErgmConstraint function should check.
    ergm_Init_abort("The CondB1Degree proposal function does not work with a non-bipartite network.")
  
  proposal
}

InitErgmProposal.CondB2Degree <- function(arguments, nw) {
  proposal <- list(name = "CondB2Degree", inputs=NULL)
  if (!is.bipartite(nw)) # Really, this should never trigger, since the InitErgmConstraint function should check.
    ergm_Init_abort("The CondB2Degree proposal function does not work with a non-bipartite network.")
  proposal
}

InitErgmProposal.CondDegreeDist <- function(arguments, nw) {
  proposal <- list(name = "CondDegreeDist", inputs=NULL)
  if (is.directed(nw)) {
    ergm_Init_warn("Using the 'degreedist' constraint with a directed network ",
          "is currently perilous.  We recommend that you use 'outdegree' or ",
          "'idegrees' instead.")
  }
  if(is.bipartite(nw)){
     proposal$name <- "BipartiteCondDegreeDist"
  }
  proposal
}

InitErgmProposal.CondInDegreeDist <- function(arguments, nw) {
  proposal <- list(name = "CondInDegreeDist", inputs=NULL)
  if (!is.directed(nw)) {
    ergm_Init_warn("Using the 'idegreedist' constraint with an undirected network ",
          "is currently perilous.  We recommend that you use 'degreedist' ",
          " instead.")
  }
  if(is.bipartite(nw)){
     proposal$name <- "BipartiteCondDegreeDist"
  }
  proposal
}

InitErgmProposal.CondOutDegreeDist <- function(arguments, nw) {
  proposal <- list(name = "CondOutDegreeDist", inputs=NULL)
  if (!is.directed(nw)) {
    ergm_Init_warn("Using the 'odegreedist' constraint with an undirected network n",
          "is currently perilous.  We recommend that you use 'degreedist' ",
          " instead.")
  }
  if(is.bipartite(nw)){
     proposal$name <- "BipartiteCondDegreeDist"
  }
  proposal
}

InitErgmProposal.ConstantEdges <- function(arguments, nw) {
  proposal <- list(name = "ConstantEdges", inputs=NULL)
  proposal
}

InitErgmProposal.HammingConstantEdges <- function(arguments, nw) {
  proposal <- list(name = "HammingConstantEdges", inputs=NULL)
  if(is.bipartite(nw)){
    proposal$name <- "BipartiteHammingConstantEdges"
  }
  proposal
}

InitErgmProposal.HammingTNT <- function(arguments, nw) {
  proposal <- list(name = "HammingTNT", inputs=NULL)
  if(is.bipartite(nw)){
    proposal$name <- "BipartiteHammingTNT"
  }
  proposal
}

InitErgmProposal.randomtoggleNonObserved <- function(arguments, nw) {
  if(network.naedgecount(nw)==0){
   ergm_Init_abort("The passed network does not have any non-observed dyads.\n Hence constraining to the observed will hold the network fixed at this network.\n Either the network or the constraint need to be altered.")
  }
  proposal <- list(name = "randomtoggleList", inputs=to_ergm_Cdouble(is.na(nw)))
  proposal
}

InitErgmProposal.NonObservedTNT <- function(arguments, nw) {
  if(network.naedgecount(nw)==0){
   ergm_Init_abort("The passed network does not have any non-observed dyads.\n Hence constraining to the observed will hold the network fixed at this network.\n Either the network or the constraint need to be altered.")
  }
  proposal <- list(name = "listTNT", inputs=to_ergm_Cdouble(is.na(nw)))
  proposal
}


InitErgmProposal.fixedas <- function(arguments, nw){
	y0<-as.edgelist(arguments$constraints$fixedas$free_dyads, prototype=nw)
	## Given the list of toggleable dyads, no formation-specific proposal function is needed:
	proposal <- list(name = "randomtoggleList", inputs=to_ergm_Cdouble(y0), pkgname="ergm")
	
	proposal
	
}

InitErgmProposal.fixedasTNT <- function(arguments, nw){
	y0<-as.edgelist(arguments$constraints$fixedas$free_dyads, prototype=nw)
	## Given the list of toggleable dyads, no formation-specific proposal function is needed:
	proposal <- list(name = "listTNT", inputs=to_ergm_Cdouble(y0), pkgname="ergm")
	
	proposal
	
}

InitErgmProposal.fixallbut <- function(arguments, nw){
	y0<-as.edgelist(arguments$constraints$fixallbut$free_dyads, prototype=nw)
	## Given the list of toggleable dyads, no formation-specific proposal function is needed:
	proposal <- list(name = "randomtoggleList", inputs=to_ergm_Cdouble(y0), pkgname="ergm")
	
	proposal
	
}


InitErgmProposal.fixallbutTNT <- function(arguments, nw){
	y0<-as.edgelist(arguments$constraints$fixallbut$free_dyads, prototype=nw)
	## Given the list of toggleable dyads, no formation-specific proposal function is needed:
	proposal <- list(name = "listTNT", inputs=to_ergm_Cdouble(y0), pkgname="ergm")
	
	proposal
	
}


InitErgmProposal.RLE <- function(arguments, nw){
  proposal <- list(name = "RLE", inputs=to_ergm_Cdouble(as.rlebdm(arguments$constraints)), pkgname="ergm")
}

InitErgmProposal.RLETNT <- function(arguments, nw){
  proposal <- list(name = "RLETNT", inputs=to_ergm_Cdouble(as.rlebdm(arguments$constraints)), pkgname="ergm")
}
