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
  if(is.directed(nw)) {
    ergm_Init_abort("BDStratTNT does not support directed networks")
  }

  # bound defaults to network.size - 1, which is effectively no bound (could be made smaller in the bipartite case, but oh well)
  bound <- NVL(arguments$bound, arguments$constraints$bd$maxout, network.size(nw) - 1L)

  # if blocks has not already been initialized, or if related arguments are passed directly to the proposal, (re)initialize it now
  if(is.null(arguments$constraints$blocks) || any(!unlist(lapply(arguments[c("blocks_attr", "levels", "levels2", "b1levels", "b2levels")], is.null)))) {
    arguments$constraints$blocks <- InitErgmConstraint.blocks(nw, attr = arguments[["blocks_attr"]], levels = arguments[["levels"]], levels2 = NVL(arguments[["levels2"]], FALSE), b1levels = arguments[["b1levels"]], b2levels = arguments[["b2levels"]])
  }

  # if Strat has not already been initialized, or if related arguments are passed directly to the proposal, (re)initialize it now
  if(is.null(arguments$constraints$Strat) || any(!unlist(lapply(arguments[c("Strat_attr", "pmat", "empirical")], is.null)))) {
    arguments$constraints$Strat <- InitErgmConstraint.Strat(nw, attr = arguments[["Strat_attr"]], pmat = arguments[["pmat"]], empirical = arguments[["empirical"]])
  }
  
  allowed.tails <- arguments$constraints$blocks$allowed.tails
  allowed.heads <- arguments$constraints$blocks$allowed.heads
  
  bd_offdiag_pairs <- which(allowed.tails != allowed.heads)  
  
  bd_tails <- c(allowed.tails, if(!is.bipartite(nw)) allowed.heads[bd_offdiag_pairs])
  bd_heads <- c(allowed.heads, if(!is.bipartite(nw)) allowed.tails[bd_offdiag_pairs])

  ## number of BD mixtypes that need to be considered when Strat mixing type is off-diag and on-diag, respectively
  bd_mixtypes <- c(length(bd_tails), length(allowed.tails))
    
  # for economy of C space, best to count # of nodes of each bd-strat pairing
  nodecountsbypairedcode <- as.integer(table(from=factor(arguments$constraints$blocks$nodecov, levels=seq_len(arguments$constraints$blocks$nlevels)), to=factor(arguments$constraints$Strat$nodecov, levels=seq_len(arguments$constraints$Strat$nlevels))))
  
  proposal <- list(name = "BDStratTNT",
                   inputs = NULL, # passed by name below
                   nmixtypes = as.integer(arguments$constraints$Strat$nmixtypes),
                   strattailattrs = as.integer(arguments$constraints$Strat$tailattrs - 1L),
                   stratheadattrs = as.integer(arguments$constraints$Strat$headattrs - 1L),
                   probvec = as.double(arguments$constraints$Strat$probvec),
                   nattrcodes = as.integer(arguments$constraints$Strat$nlevels),
                   strat_vattr = as.integer(arguments$constraints$Strat$nodecov - 1L),
                   indmat = as.integer(t(arguments$constraints$Strat$indmat)), 
                   nodecountsbypairedcode = as.integer(nodecountsbypairedcode),
                   bound = as.integer(bound),
                   bd_levels = as.integer(arguments$constraints$blocks$nlevels),
                   bd_vattr = as.integer(arguments$constraints$blocks$nodecov - 1L),
                   bd_tails = as.integer(bd_tails - 1L),
                   bd_heads = as.integer(bd_heads - 1L),
                   bd_mixtypes = as.integer(bd_mixtypes),
                   empirical_flag = as.integer(arguments$constraints$Strat$empirical),
                   amat_C = as.integer(t(arguments$constraints$blocks$amat_C)))

  proposal
}

InitErgmProposal.BDTNT <- function(arguments, nw) {
  # BDTNT does not currently support directed networks
  if(is.directed(nw)) {
    ergm_Init_abort(sQuote("BDTNT"), " only supports undirected networks.")
  }
  
  # bound defaults to network.size - 1, which is effectively no bound (could be made smaller in the bipartite case, but oh well)
  bound <- NVL(arguments$bound, arguments$constraints$bd$maxout, network.size(nw) - 1L)

  # if blocks has not already been initialized, or if related arguments are passed directly to the proposal, (re)initialize it now
  if(is.null(arguments$constraints$blocks) || any(!unlist(lapply(arguments[c("attr", "levels", "levels2", "b1levels", "b2levels")], is.null)))) {
    arguments$constraints$blocks <- InitErgmConstraint.blocks(nw, attr = arguments[["attr"]], levels = arguments[["levels"]], levels2 = NVL(arguments[["levels2"]], FALSE), b1levels = arguments[["b1levels"]], b2levels = arguments[["b2levels"]])
  }

  proposal <- list(name = "BDTNT", 
                   inputs = NULL, # passed by name below
                   bound = as.integer(bound),
                   nlevels = as.integer(arguments$constraints$blocks$nlevels),
                   nodecountsbycode = as.integer(arguments$constraints$blocks$nodecountsbycode),
                   nmixtypes = as.integer(arguments$constraints$blocks$nmixtypes),
                   allowed.tails = as.integer(arguments$constraints$blocks$allowed.tails - 1L),
                   allowed.heads = as.integer(arguments$constraints$blocks$allowed.heads - 1L),
                   nodecov = as.integer(arguments$constraints$blocks$nodecov - 1L),
                   amat_C = as.integer(t(arguments$constraints$blocks$amat_C)))
                   
  proposal
}

InitErgmProposal.StratTNT <- function(arguments, nw) {
  # if Strat has not already been initialized, or if related arguments are passed directly to the proposal, (re)initialize it now
  if(is.null(arguments$constraints$Strat) || any(!unlist(lapply(arguments[c("attr", "pmat", "empirical")], is.null)))) {
    arguments$constraints$Strat <- InitErgmConstraint.Strat(nw, attr = arguments[["attr"]], pmat = arguments[["pmat"]], empirical = arguments[["empirical"]])
  }

  ## subtract one from attr codes for greater convenience re. C's zero-based indexing  
  proposal <- list(name = "StratTNT", 
                   inputs = NULL, # passed by name below
                   nmixtypes = as.integer(arguments$constraints$Strat$nmixtypes),
                   tailattrs = as.integer(arguments$constraints$Strat$tailattrs - 1L),
                   headattrs = as.integer(arguments$constraints$Strat$headattrs - 1L),
                   probvec = as.double(arguments$constraints$Strat$probvec),
                   nlevels = as.integer(arguments$constraints$Strat$nlevels),
                   nodecountsbycode = as.integer(arguments$constraints$Strat$nodecountsbycode),
                   nodeindicesbycode = as.integer(arguments$constraints$Strat$nodeindicesbycode),
                   nodecov = as.integer(arguments$constraints$Strat$nodecov - 1L),
                   indmat = as.integer(t(arguments$constraints$Strat$indmat)),
                   empirical = as.integer(arguments$constraints$Strat$empirical))

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

InitErgmProposal.NonObservedTNT <- function(arguments, nw) {
  if(network.naedgecount(nw)==0){
   ergm_Init_abort("The passed network does not have any non-observed dyads.\n Hence constraining to the observed will hold the network fixed at this network.\n Either the network or the constraint need to be altered.")
  }
  proposal <- list(name = "listTNT", iinputs=to_ergm_Cdouble(is.na(nw)))
  proposal
}

InitErgmProposal.fixedasTNT <- function(arguments, nw){
	y0<-as.edgelist(arguments$constraints$fixedas$free_dyads, prototype=nw)
	## Given the list of toggleable dyads, no formation-specific proposal function is needed:
	proposal <- list(name = "listTNT", iinputs=to_ergm_Cdouble(y0), pkgname="ergm")
	
	proposal
	
}

InitErgmProposal.fixallbutTNT <- function(arguments, nw){
	y0<-as.edgelist(arguments$constraints$fixallbut$free_dyads, prototype=nw)
	## Given the list of toggleable dyads, no formation-specific proposal function is needed:
	proposal <- list(name = "listTNT", iinputs=to_ergm_Cdouble(y0), pkgname="ergm")
	
	proposal
	
}

InitErgmProposal.RLETNT <- function(arguments, nw){
  proposal <- list(name = "RLETNT", inputs=to_ergm_Cdouble(as.rlebdm(arguments$constraints)), pkgname="ergm")
}
