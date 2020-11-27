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
  arguments$bound <- NVL(arguments$bound, arguments$constraints$BD$bound)
  arguments$BD_attr <- NVL(arguments$BD_attr, arguments$constraints$BD$attr)
  arguments$fmat <- NVL(arguments$fmat, arguments$constraints$BD$fmat)

  arguments$Strat_attr <- NVL(arguments$Strat_attr, arguments$constraints$Strat$attr)
  arguments$pmat <- NVL(arguments$pmat, arguments$constraints$Strat$pmat)
  arguments$empirical <- NVL(arguments$empirical, arguments$constraints$Strat$empirical)

  if(is.directed(nw)) {
    ergm_Init_abort("BDStratTNT does not support directed networks")
  }
  
  if(is.bipartite(nw)) {
    strat_row_nodecov <- NVL2(arguments$Strat_attr, ergm_get_vattr(arguments$Strat_attr, nw, bip="b1"), rep(1, nw %n% "bipartite"))
    strat_col_nodecov <- NVL2(arguments$Strat_attr, ergm_get_vattr(arguments$Strat_attr, nw, bip="b2"), rep(1, network.size(nw) - (nw %n% "bipartite")))

    strat_row_levels <- sort(unique(strat_row_nodecov))
    strat_col_levels <- sort(unique(strat_col_nodecov))

    strat_levels <- c(strat_row_levels, strat_col_levels)

    strat_row_nodecov <- match(strat_row_nodecov, strat_row_levels)
    strat_col_nodecov <- match(strat_col_nodecov, strat_col_levels)
    strat_nodecov <- c(strat_row_nodecov, strat_col_nodecov + length(strat_row_levels))
  } else {
    strat_row_nodecov <- NVL2(arguments$Strat_attr, ergm_get_vattr(arguments$Strat_attr, nw), rep(1, network.size(nw)))
    strat_col_nodecov <- strat_row_nodecov

    strat_row_levels <- sort(unique(strat_row_nodecov))
    strat_col_levels <- strat_row_levels

    strat_levels <- strat_row_levels
    
    strat_row_nodecov <- match(strat_row_nodecov, strat_row_levels)
    strat_col_nodecov <- strat_row_nodecov
    strat_nodecov <- strat_row_nodecov
  }
  
  pmat <- NVL(arguments$pmat, matrix(1, nrow = length(strat_row_levels), ncol = length(strat_col_levels)))
    
  if(NROW(pmat) != length(strat_row_levels) || NCOL(pmat) != length(strat_col_levels)) {
    ergm_Init_abort(sQuote("pmat"), " does not have the correct dimensions for ", sQuote("Strat_attr"), ".")    
  }
  
  if(!is.bipartite(nw)) {
    # for undirected unipartite, symmetrize pmat and then set the sub-diagonal to zero
    pmat <- (pmat + t(pmat))/2
    pmat[lower.tri(pmat)] <- 0
  }  
    
  # renormalize to probability matrix
  pmat <- pmat/sum(pmat)

  # record the tail and head attr code for each mixing type with positive probability
  prob_inds <- which(pmat > 0, arr.ind = TRUE)
  tailattrs <- prob_inds[,1]
  headattrs <- prob_inds[,2]
  probvec <- pmat[prob_inds]
  
  if(is.bipartite(nw)) {
    headattrs <- headattrs + length(strat_row_levels)
  }
      
  bound <- NVL(arguments$bound, network.size(nw) - 1)
  
  if(is.bipartite(nw)) {
    bd_row_nodecov <- NVL2(arguments$BD_attr, ergm_get_vattr(arguments$BD_attr, nw, bip="b1"), rep(1, nw %n% "bipartite"))
    bd_col_nodecov <- NVL2(arguments$BD_attr, ergm_get_vattr(arguments$BD_attr, nw, bip="b2"), rep(1, network.size(nw) - (nw %n% "bipartite")))

    bd_row_levels <- sort(unique(bd_row_nodecov))
    bd_col_levels <- sort(unique(bd_col_nodecov))

    bd_levels <- c(bd_row_levels, bd_col_levels)

    bd_row_nodecov <- match(bd_row_nodecov, bd_row_levels)
    bd_col_nodecov <- match(bd_col_nodecov, bd_col_levels)
    bd_nodecov <- c(bd_row_nodecov, bd_col_nodecov + length(bd_row_levels))
  } else {
    bd_row_nodecov <- NVL2(arguments$BD_attr, ergm_get_vattr(arguments$BD_attr, nw), rep(1, network.size(nw)))
    bd_col_nodecov <- bd_row_nodecov

    bd_row_levels <- sort(unique(bd_row_nodecov))
    bd_col_levels <- bd_row_levels

    bd_levels <- bd_row_levels
    
    bd_row_nodecov <- match(bd_row_nodecov, bd_row_levels)
    bd_col_nodecov <- bd_row_nodecov
    bd_nodecov <- bd_row_nodecov
  }
  
  # by default, no pairings are forbidden
  fmat <- NVL(arguments$fmat, matrix(FALSE, nrow = length(bd_row_levels), ncol = length(bd_col_levels)))
  
  if(NROW(fmat) != length(bd_row_levels) || NCOL(fmat) != length(bd_col_levels)) {
    ergm_Init_abort(sQuote("fmat"), " does not have the correct dimensions for ", sQuote("BD_attr"), ".")
  }    
  
  if(!is.bipartite(nw)) {
    # for undirected unipartite, symmetrize fmat and then set the sub-diagonal to TRUE
    fmat <- fmat | t(fmat)
    fmat[lower.tri(fmat)] <- TRUE
  }
  
  # create vectors of allowed mixing types
  allowed.attrs <- which(!fmat, arr.ind = TRUE)
  allowed.tails <- allowed.attrs[,1]
  allowed.heads <- allowed.attrs[,2]
  
  if(is.bipartite(nw)) {
    allowed.heads <- allowed.heads + length(bd_row_levels)  
  }
 
  bd_offdiag_pairs <- which(allowed.tails != allowed.heads)
  
  type1_bd_tails <- allowed.tails
  type1_bd_heads <- allowed.heads
  
  type2_bd_tails <- c(allowed.tails, allowed.heads[bd_offdiag_pairs])
  type2_bd_heads <- c(allowed.heads, allowed.tails[bd_offdiag_pairs])

  BDtailsbyStrattype <- vector(mode = "list", length = length(tailattrs))
  BDheadsbyStrattype <- vector(mode = "list", length = length(tailattrs))
  
  for(i in 1:length(tailattrs)) {
    if(tailattrs[i] == headattrs[i] || is.bipartite(nw)) {
      BDtailsbyStrattype[[i]] <- type1_bd_tails
      BDheadsbyStrattype[[i]] <- type1_bd_heads
    } else {
      BDtailsbyStrattype[[i]] <- type2_bd_tails
      BDheadsbyStrattype[[i]] <- type2_bd_heads    
    }
  }

  BDtypesbyStrattype <- sapply(BDtailsbyStrattype, length)
  BDtailsbyStrattype <- unlist(BDtailsbyStrattype)
  BDheadsbyStrattype <- unlist(BDheadsbyStrattype)
  
  indmat <- matrix(-1L, nrow=length(strat_levels), ncol=length(strat_levels))
  indmat[cbind(tailattrs, headattrs)] <- seq_along(tailattrs) - 1L  # zero-based for C code
  if(!is.bipartite(nw)) {
    # symmetrize for undirected unipartite
    indmat[cbind(headattrs, tailattrs)] <- seq_along(tailattrs) - 1L
  }

  # for economy of C space, best to count # of nodes of each bd-strat pairing
  nodecountsbypairedcode <- as.integer(table(from=bd_nodecov, to=strat_nodecov))
    
  empirical_flag <- as.logical(NVL(arguments$empirical, FALSE))

  proposal <- list(name = "BDStratTNT",
                   inputs = NULL, # passed by name below
                   nmixtypes = as.integer(length(tailattrs)),
                   strattailattrs = as.integer(tailattrs - 1L),
                   stratheadattrs = as.integer(headattrs - 1L),
                   probvec = as.double(probvec),
                   nattrcodes = as.integer(length(strat_levels)),
                   strat_vattr = as.integer(strat_nodecov - 1L),
                   indmat = as.integer(t(indmat)), 
                   nodecountsbypairedcode = as.integer(nodecountsbypairedcode),
                   bound = as.integer(bound),
                   bd_levels = as.integer(length(bd_levels)),
                   bd_vattr = as.integer(bd_nodecov - 1L),
                   BDtypesbyStrattype = as.integer(BDtypesbyStrattype), 
                   BDtailsbyStrattype = as.integer(BDtailsbyStrattype - 1L), 
                   BDheadsbyStrattype = as.integer(BDheadsbyStrattype - 1L), 
                   empirical_flag = as.integer(empirical_flag))

  proposal
}

InitErgmProposal.BDTNT <- function(arguments, nw) {
  arguments$bound <- NVL(arguments$bound, arguments$constraints$BD$bound)
  arguments$attr <- NVL(arguments$attr, arguments$constraints$BD$attr)
  arguments$fmat <- NVL(arguments$fmat, arguments$constraints$BD$fmat)

  # BDTNT does not currently support directed networks
  if(is.directed(nw)) {
    ergm_Init_abort(sQuote("BDTNT"), " only supports undirected networks.")
  }
  
  if(is.bipartite(nw)) {
    bd_row_nodecov <- NVL2(arguments$attr, ergm_get_vattr(arguments$attr, nw, bip="b1"), rep(1, nw %n% "bipartite"))
    bd_col_nodecov <- NVL2(arguments$attr, ergm_get_vattr(arguments$attr, nw, bip="b2"), rep(1, network.size(nw) - (nw %n% "bipartite")))

    bd_row_levels <- sort(unique(bd_row_nodecov))
    bd_col_levels <- sort(unique(bd_col_nodecov))

    bd_levels <- c(bd_row_levels, bd_col_levels)

    bd_row_nodecov <- match(bd_row_nodecov, bd_row_levels)
    bd_col_nodecov <- match(bd_col_nodecov, bd_col_levels)
    bd_nodecov <- c(bd_row_nodecov, bd_col_nodecov + length(bd_row_levels))
  } else {
    bd_row_nodecov <- NVL2(arguments$attr, ergm_get_vattr(arguments$attr, nw), rep(1, network.size(nw)))
    bd_col_nodecov <- bd_row_nodecov

    bd_row_levels <- sort(unique(bd_row_nodecov))
    bd_col_levels <- bd_row_levels

    bd_levels <- bd_row_levels
    
    bd_row_nodecov <- match(bd_row_nodecov, bd_row_levels)
    bd_col_nodecov <- bd_row_nodecov
    bd_nodecov <- bd_row_nodecov
  }
  
  # by default, no pairings are forbidden
  fmat <- NVL(arguments$fmat, matrix(FALSE, nrow = length(bd_row_levels), ncol = length(bd_col_levels)))
  
  if(NROW(fmat) != length(bd_row_levels) || NCOL(fmat) != length(bd_col_levels)) {
    ergm_Init_abort(sQuote("fmat"), " does not have the correct dimensions for ", sQuote("BD_attr"), ".")
  }    
  
  if(!is.bipartite(nw)) {
    # for undirected unipartite, symmetrize fmat and then set the sub-diagonal to TRUE
    fmat <- fmat | t(fmat)
    fmat[lower.tri(fmat)] <- TRUE
  }
  
  # create vectors of allowed mixing types
  allowed.attrs <- which(!fmat, arr.ind = TRUE)
  allowed.tails <- allowed.attrs[,1]
  allowed.heads <- allowed.attrs[,2]
  
  if(is.bipartite(nw)) {
    allowed.heads <- allowed.heads + length(bd_row_levels)  
  }
  
  # bound defaults to network.size - 1, which is effectively no bound (could be made smaller in the bipartite case, but oh well)
  bound <- NVL(arguments$bound, network.size(nw) - 1L)  
  
  ncodes <- length(bd_levels)
  
  # record number of nodes of each type
  nodecountsbycode <- tabulate(bd_nodecov, nbins=length(bd_levels))
  
  ## subtract one from attr codes for greater convenience re. C's zero-based indexing  
  proposal <- list(name = "BDTNT", 
                   inputs = NULL, # passed by name below
                   bound = as.integer(bound),
                   nlevels = as.integer(ncodes),
                   nodecountsbycode = as.integer(nodecountsbycode),
                   nmixtypes = length(allowed.tails),
                   allowed.tails = as.integer(allowed.tails - 1L),
                   allowed.heads = as.integer(allowed.heads - 1L),
                   nodecov = as.integer(bd_nodecov - 1L))
                   
  proposal
}

InitErgmProposal.StratTNT <- function(arguments, nw) {
  arguments$attr <- NVL(arguments$attr, arguments$constraints$Strat$attr)
  arguments$pmat <- NVL(arguments$pmat, arguments$constraints$Strat$pmat)
  arguments$empirical <- NVL(arguments$empirical, arguments$constraints$Strat$empirical)
  
  if(is.bipartite(nw)) {
    strat_row_nodecov <- NVL2(arguments$attr, ergm_get_vattr(arguments$attr, nw, bip="b1"), rep(1, nw %n% "bipartite"))
    strat_col_nodecov <- NVL2(arguments$attr, ergm_get_vattr(arguments$attr, nw, bip="b2"), rep(1, network.size(nw) - (nw %n% "bipartite")))

    strat_row_levels <- sort(unique(strat_row_nodecov))
    strat_col_levels <- sort(unique(strat_col_nodecov))

    strat_levels <- c(strat_row_levels, strat_col_levels)

    strat_row_nodecov <- match(strat_row_nodecov, strat_row_levels)
    strat_col_nodecov <- match(strat_col_nodecov, strat_col_levels)
    strat_nodecov <- c(strat_row_nodecov, strat_col_nodecov + length(strat_row_levels))
  } else {
    strat_row_nodecov <- NVL2(arguments$attr, ergm_get_vattr(arguments$attr, nw), rep(1, network.size(nw)))
    strat_col_nodecov <- strat_row_nodecov

    strat_row_levels <- sort(unique(strat_row_nodecov))
    strat_col_levels <- strat_row_levels

    strat_levels <- strat_row_levels
    
    strat_row_nodecov <- match(strat_row_nodecov, strat_row_levels)
    strat_col_nodecov <- strat_row_nodecov
    strat_nodecov <- strat_row_nodecov
  }
  
  pmat <- NVL(arguments$pmat, matrix(1, nrow = length(strat_row_levels), ncol = length(strat_col_levels)))
    
  if(NROW(pmat) != length(strat_row_levels) || NCOL(pmat) != length(strat_col_levels)) {
    ergm_Init_abort(sQuote("pmat"), " does not have the correct dimensions for ", sQuote("Strat_attr"), ".")    
  }
  
  if(!is.bipartite(nw) && !is.directed(nw)) {
    # for undirected unipartite, symmetrize pmat and then set the sub-diagonal to zero
    pmat <- (pmat + t(pmat))/2
    pmat[lower.tri(pmat)] <- 0
  }  
    
  # renormalize to probability matrix
  pmat <- pmat/sum(pmat)

  # record the tail and head attr code for each mixing type with positive probability
  prob_inds <- which(pmat > 0, arr.ind = TRUE)
  tailattrs <- prob_inds[,1]
  headattrs <- prob_inds[,2]
  probvec <- pmat[prob_inds]
  
  if(is.bipartite(nw)) {
    headattrs <- headattrs + length(strat_row_levels)
  }

  probvec <- cumsum(probvec)

  # record the number of mixing types and the number of unique attr codes
  nmixingtypes <- length(probvec)
  ncodes <- length(strat_levels)
    
  # record the nodal indices grouped (and counted) by attr code
  nodeindicesbycode <- order(strat_nodecov)
  nodecountsbycode <- tabulate(strat_nodecov, nbins=length(strat_levels))
  
  ## may wish to add check that if pmat[i,j] > 0 then at least one dyad  
  indmat <- matrix(-1L, nrow=length(strat_levels), ncol=length(strat_levels))
  indmat[cbind(tailattrs, headattrs)] <- seq_along(tailattrs) - 1L  # zero-based for C code
  if(!is.bipartite(nw) && !is.directed(nw)) {
    # symmetrize for undirected unipartite
    indmat[cbind(headattrs, tailattrs)] <- seq_along(tailattrs) - 1L
  }
  
  empirical_flag <- as.logical(NVL(arguments$empirical, FALSE))

  ## subtract one from attr codes for greater convenience re. C's zero-based indexing  
  proposal <- list(name = "StratTNT", 
                   inputs = NULL, # passed by name below
                   nmixtypes = as.integer(nmixingtypes),
                   tailattrs = as.integer(tailattrs - 1L),
                   headattrs = as.integer(headattrs - 1L),
                   probvec = as.double(probvec),
                   nlevels = as.integer(ncodes),
                   nodecountsbycode = as.integer(nodecountsbycode),
                   nodeindicesbycode = as.integer(nodeindicesbycode),
                   nodecov = as.integer(strat_nodecov - 1L),
                   indmat = as.integer(t(indmat)),
                   empirical = as.integer(empirical_flag))

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
