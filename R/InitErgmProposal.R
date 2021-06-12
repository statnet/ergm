#  File R/InitErgmProposal.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################
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
DyadGenType <- list(RandDyadGen=0L, WtRandDyadGen=1L, RLEBDM1DGen=2L, EdgeListGen=3L)

InitErgmProposal.randomtoggle <- function(arguments, nw){
  list(name = "randomtoggle", dyadgen = ergm_dyadgen_select(arguments, nw), bd = ergm_bd_init(arguments, nw))
}

InitErgmProposal.TNT <- function(nw, arguments, ...){
  list(name = "TNT", dyadgen = ergm_dyadgen_select(arguments, nw), bd = ergm_bd_init(arguments, nw))
}

InitErgmProposal.BDStratTNT <- function(arguments, nw) {
  # if bd has not already been initialized, or if related arguments are passed directly to the proposal, (re)initialize it now
  if(is.null(arguments$constraints$bd) || any(!unlist(lapply(arguments[c("attribs", "maxout", "maxin")], is.null)))) {
    arguments$constraints$bd <- InitErgmConstraint.bd(nw, list(attribs = arguments[["attribs"]], maxout = arguments[["maxout"]], maxin = arguments[["maxin"]]))
  }

  attribs <- NVL(arguments$constraints$bd$attribs, matrix(TRUE, ncol = 1L, nrow = network.size(nw)))
  
  maxout <- NVL(arguments$constraints$bd$maxout, network.size(nw))
  maxout[is.na(maxout)] <- network.size(nw)
  maxout <- matrix(rep(maxout, length.out = length(attribs)), ncol = ncol(attribs))
  
  maxin <- NVL(arguments$constraints$bd$maxin, network.size(nw))
  maxin[is.na(maxin)] <- network.size(nw)
  maxin <- matrix(rep(maxin, length.out = length(attribs)), ncol = ncol(attribs))

  bd_vattr <- which(attribs, arr.ind = TRUE)
  bd_vattr <- bd_vattr[order(bd_vattr[,1L]), 2L]
    
  ## which attribute pairings are allowed by bd? only consider maxin if directed
  bd_mixmat <- matrix(FALSE, nrow = NCOL(attribs), ncol = NCOL(attribs))
  
  maxout_pairs <- which(maxout > 0, arr.ind = TRUE)
  maxout_pairs[,1L] <- bd_vattr[maxout_pairs[,1L]]
  bd_mixmat[maxout_pairs] <- TRUE
  if(is.directed(nw)) {
    maxin_pairs <- which(maxin > 0, arr.ind = TRUE)
    maxin_pairs[,1L] <- bd_vattr[maxin_pairs[,1L]]
    bd_mixmat[maxin_pairs[,c(2L,1L)]] <- TRUE
  } else {
    bd_mixmat <- bd_mixmat | t(bd_mixmat)
    bd_mixmat[lower.tri(bd_mixmat, diag = FALSE)] <- FALSE
  }

  bd_pairs <- which(bd_mixmat, arr.ind = TRUE)
  bd_tails <- bd_pairs[,1L]
  bd_heads <- bd_pairs[,2L]
  bd_nlevels <- NCOL(attribs)

  ## need to handle undirected case as for blocks, but bipartite along with unip for now
  bd_offdiag_pairs <- which(bd_tails != bd_heads)  
  
  bd_allowed_tails <- c(bd_tails, if(!is.directed(nw)) bd_heads[bd_offdiag_pairs])
  bd_allowed_heads <- c(bd_heads, if(!is.directed(nw)) bd_tails[bd_offdiag_pairs])

  ## number of bd mixtypes that need to be considered when strat and blocks mixing types are off-diag and on-diag, respectively
  bd_nmixtypes <- c(length(bd_allowed_tails), length(bd_tails))
  
  
  
  # if blocks has not already been initialized, or if related arguments are passed directly to the proposal, (re)initialize it now
  if(is.null(arguments$constraints$blocks) || any(!unlist(lapply(arguments[c("blocks_attr", "levels", "levels2", "b1levels", "b2levels")], is.null)))) {
    arguments$constraints$blocks <- InitErgmConstraint.blocks(nw, list(attr = arguments[["blocks_attr"]], levels = arguments[["levels"]], levels2 = NVL(arguments[["levels2"]], FALSE), b1levels = arguments[["b1levels"]], b2levels = arguments[["b2levels"]]))
  }

  # check for old name
  if(is.null(arguments$constraints$strat) && !is.null(arguments$constraints$Strat)) {
    arguments$constraints$strat <- arguments$constraints$Strat
  }

  # if strat has not already been initialized, or if related arguments are passed directly to the proposal, (re)initialize it now
  if(is.null(arguments$constraints$strat) || any(!unlist(lapply(arguments[c("strat_attr", "pmat", "empirical")], is.null)))) {
    arguments$constraints$strat <- InitErgmConstraint.strat(nw, list(attr = arguments[["strat_attr"]], pmat = arguments[["pmat"]], empirical = NVL(arguments[["empirical"]],FALSE)))
  }

  nodecov <- arguments$constraints$blocks$nodecov
  amat <- arguments$constraints$blocks$amat

  if(is.bipartite(nw)) {
    nodecov[-seq_len(nw %n% "bipartite")] <- nodecov[-seq_len(nw %n% "bipartite")] + NROW(amat)
    pairs_mat <- matrix(FALSE, nrow = NROW(amat) + NCOL(amat), ncol = NROW(amat) + NCOL(amat))
    pairs_mat[seq_len(NROW(amat)), -seq_len(NROW(amat))] <- amat
  } else if(!is.directed(nw)) {
    pairs_mat <- amat
    pairs_mat[lower.tri(pairs_mat, diag = FALSE)] <- FALSE
  } else {
    pairs_mat <- amat
  }
  
  allowed.attrs <- which(pairs_mat, arr.ind = TRUE)
  allowed.tails <- allowed.attrs[,1]
  allowed.heads <- allowed.attrs[,2]  

  nlevels <- NROW(pairs_mat)
  nodecountsbycode <- tabulate(nodecov, nbins = nlevels)
  
  if(!is.directed(nw) && !is.bipartite(nw)) {
    pairs_mat <- pairs_mat | t(pairs_mat)
  }
  
  pairs_to_keep <- (allowed.tails != allowed.heads & nodecountsbycode[allowed.tails] > 0 & nodecountsbycode[allowed.heads] > 0) | (allowed.tails == allowed.heads & nodecountsbycode[allowed.tails] > 1)
  allowed.tails <- allowed.tails[pairs_to_keep]
  allowed.heads <- allowed.heads[pairs_to_keep]
  
  
  blocks_offdiag_pairs <- which(allowed.tails != allowed.heads)  
  
  blocks_tails <- c(allowed.tails, if(!is.bipartite(nw) && !is.directed(nw)) allowed.heads[blocks_offdiag_pairs])
  blocks_heads <- c(allowed.heads, if(!is.bipartite(nw) && !is.directed(nw)) allowed.tails[blocks_offdiag_pairs])

  ## number of blocks mixtypes that need to be considered when strat mixing type is off-diag and on-diag, respectively
  blocks_mixtypes <- c(length(blocks_tails), length(allowed.tails))
    
  # for economy of C space, best to count # of nodes of each bd-strat pairing
  nodecountsbyjointcode <- as.integer(table(factor(bd_vattr, levels=seq_len(bd_nlevels)),
                                            factor(nodecov, levels=seq_len(nlevels)), 
                                            factor(arguments$constraints$strat$nodecov, levels=seq_len(arguments$constraints$strat$nlevels))))
  
  proposal <- list(name = "BDStratTNT",
                   inputs = NULL, # passed by name below
                   nmixtypes = as.integer(arguments$constraints$strat$nmixtypes),
                   strattailattrs = as.integer(arguments$constraints$strat$tailattrs - 1L),
                   stratheadattrs = as.integer(arguments$constraints$strat$headattrs - 1L),
                   probvec = as.double(arguments$constraints$strat$probvec),
                   nattrcodes = as.integer(arguments$constraints$strat$nlevels),
                   strat_vattr = as.integer(arguments$constraints$strat$nodecov - 1L),
                   indmat = as.integer(t(arguments$constraints$strat$indmat)),
                   nodecountsbyjointcode = as.integer(nodecountsbyjointcode),
                   maxout = as.integer(maxout),
                   maxin = as.integer(maxin),
                   bd_vattr = as.integer(bd_vattr - 1L),
                   bd_tails = as.integer(bd_allowed_tails - 1L),
                   bd_heads = as.integer(bd_allowed_heads - 1L),
                   bd_nlevels = as.integer(bd_nlevels),
                   bd_nmixtypes = as.integer(bd_nmixtypes),
                   blocks_levels = as.integer(nlevels),
                   blocks_vattr = as.integer(nodecov - 1L),
                   blocks_tails = as.integer(blocks_tails - 1L),
                   blocks_heads = as.integer(blocks_heads - 1L),
                   blocks_mixtypes = as.integer(blocks_mixtypes),
                   empirical_flag = as.integer(arguments$constraints$strat$empirical),
                   amat = as.integer(t(pairs_mat)))

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
  proposals$bd <- ergm_bd_init(arguments, nw)
  proposal
}

InitErgmProposal.ConstantEdges <- function(arguments, nw) {
  proposal <- list(name = "ConstantEdges", bd = ergm_bd_init(arguments, nw))
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
