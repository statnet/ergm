#  File R/InitErgmProposal.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
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
DyadGenType <- list(RandDyadGen=0L, WtRandDyadGen=1L, RLEBDM1DGen=2L, EdgeListGen=3L)

InitErgmProposal.randomtoggle <- function(arguments, nw){
  list(name = "randomtoggle", dyadgen = ergm_dyadgen_select(arguments, nw))
}

InitErgmProposal.TNT <- function(arguments, nw){
  list(name = "TNT", dyadgen = ergm_dyadgen_select(arguments, nw))
}

InitErgmProposal.BDStratTNT <- function(arguments, nw) {
  # maxout defaults to network.size - 1, which is effectively no bound (could be made smaller in the bipartite case, but oh well)
  maxout <- NVL(arguments$maxout, arguments$constraints$bd$maxout, network.size(nw) - 1L)
  if(is.na(maxout)) maxout <- network.size(nw) - 1L

  # maxin defaults to network.size - 1, which is effectively no bound (could be made smaller in the bipartite case, but oh well)
  maxin <- NVL(arguments$maxin, arguments$constraints$bd$maxin, network.size(nw) - 1L)
  if(is.na(maxin)) maxin <- network.size(nw) - 1L

  # if blocks has not already been initialized, or if related arguments are passed directly to the proposal, (re)initialize it now
  if(is.null(arguments$constraints$blocks) || any(!unlist(lapply(arguments[c("blocks_attr", "levels", "levels2", "b1levels", "b2levels")], is.null)))) {
    arguments$constraints$blocks <- InitErgmConstraint.blocks(nw, attr = arguments[["blocks_attr"]], levels = arguments[["levels"]], levels2 = NVL(arguments[["levels2"]], FALSE), b1levels = arguments[["b1levels"]], b2levels = arguments[["b2levels"]])
  }

  # check for old name
  if(is.null(arguments$constraints$strat) && !is.null(arguments$constraints$Strat)) {
    arguments$constraints$strat <- arguments$constraints$Strat
  }

  # if strat has not already been initialized, or if related arguments are passed directly to the proposal, (re)initialize it now
  if(is.null(arguments$constraints$strat) || any(!unlist(lapply(arguments[c("strat_attr", "pmat", "empirical")], is.null)))) {
    arguments$constraints$strat <- InitErgmConstraint.strat(nw, attr = arguments[["strat_attr"]], pmat = arguments[["pmat"]], empirical = arguments[["empirical"]])
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
  
  
  bd_offdiag_pairs <- which(allowed.tails != allowed.heads)  
  
  bd_tails <- c(allowed.tails, if(!is.bipartite(nw) && !is.directed(nw)) allowed.heads[bd_offdiag_pairs])
  bd_heads <- c(allowed.heads, if(!is.bipartite(nw) && !is.directed(nw)) allowed.tails[bd_offdiag_pairs])

  ## number of BD mixtypes that need to be considered when strat mixing type is off-diag and on-diag, respectively
  bd_mixtypes <- c(length(bd_tails), length(allowed.tails))
    
  # for economy of C space, best to count # of nodes of each bd-strat pairing
  nodecountsbypairedcode <- as.integer(table(from=factor(nodecov, levels=seq_len(nlevels)), to=factor(arguments$constraints$strat$nodecov, levels=seq_len(arguments$constraints$strat$nlevels))))
  
  proposal <- list(name = "BDStratTNT",
                   inputs = NULL, # passed by name below
                   nmixtypes = as.integer(arguments$constraints$strat$nmixtypes),
                   strattailattrs = as.integer(arguments$constraints$strat$tailattrs - 1L),
                   stratheadattrs = as.integer(arguments$constraints$strat$headattrs - 1L),
                   probvec = as.double(arguments$constraints$strat$probvec),
                   nattrcodes = as.integer(arguments$constraints$strat$nlevels),
                   strat_vattr = as.integer(arguments$constraints$strat$nodecov - 1L),
                   indmat = as.integer(t(arguments$constraints$strat$indmat)),
                   nodecountsbypairedcode = as.integer(nodecountsbypairedcode),
                   maxout = as.integer(maxout),
                   maxin = as.integer(maxin),
                   bd_levels = as.integer(nlevels),
                   bd_vattr = as.integer(nodecov - 1L),
                   bd_tails = as.integer(bd_tails - 1L),
                   bd_heads = as.integer(bd_heads - 1L),
                   bd_mixtypes = as.integer(bd_mixtypes),
                   empirical_flag = as.integer(arguments$constraints$strat$empirical),
                   amat = as.integer(t(pairs_mat)),
                   skip_bd = TRUE)

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
