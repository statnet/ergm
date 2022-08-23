#  File R/InitErgmProposal.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
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

#' @templateVar name randomtoggle
#' @aliases InitErgmProposal.randomtoggle
#' @title Propose a randomly selected dyad to toggle
#' @description Propose a randomly selected dyad to toggle
#' @template ergmProposal-general
NULL
InitErgmProposal.randomtoggle <- function(arguments, nw){
  list(name = "randomtoggle", dyadgen = ergm_dyadgen_select(arguments, nw), bd = ergm_bd_init(arguments, nw))
}

#' @templateVar name TNT
#' @aliases InitErgmProposal.TNT
#' @title Default MH algorithm
#' @description Stratifies the population of dyads
#'   edge status: those having ties and those having no ties (hence T/NT).
#'   This is useful for improving performance in sparse networks,
#'   because it gives at least 50\% chance of proposing a toggle of an existing edge.
#' @template ergmProposal-general
NULL
InitErgmProposal.TNT <- function(nw, arguments, ...){
  list(name = "TNT", dyadgen = ergm_dyadgen_select(arguments, nw), bd = ergm_bd_init(arguments, nw))
}

#' @templateVar name BDStratTNT
#' @aliases InitErgmProposal.BDStratTNT
#' @title TNT proposal with degree bounds
#' @description Implements a TNT proposal with any subset of the following features:
#'   1. upper on degree, specified via the [`bd`][bd-ergmConstraint] constraint's `maxout`, `maxin`, and `attribs` arguments
#'   2. stratification of proposals according to mixing type on a vertex attribute, specified via the [`strat`][strat-ergmHint] hint;
#'   3. fixation of specified mixing types on a(nother) vertex attribute, specified via the [`blocks`][blocks-ergmConstraint] constraint.
#' @template ergmProposal-general
NULL
InitErgmProposal.BDStratTNT <- function(arguments, nw) {
  # if bd has not already been initialized, or if related arguments are passed directly to the proposal, (re)initialize it now
  if(any(!unlist(lapply(arguments[c("attribs", "maxout", "maxin")], is.null)))) {
    arguments$constraints$bd <- InitErgmConstraint.bd(nw, list(attribs = arguments[["attribs"]], maxout = arguments[["maxout"]], maxin = arguments[["maxin"]]))
  }

  dyad_indep <- is.null(arguments$constraints$bd)
  
  if(is.null(arguments$constraints$bd)) {
    arguments$constraints$bd <- InitErgmConstraint.bd(nw, list())
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
  if(length(intersect(names(arguments), c("blocks_attr", "levels", "levels2", "b1levels", "b2levels"))) > 0) {
    arglist <- arguments[intersect(names(arguments), c("blocks_attr", "levels", "levels2", "b1levels", "b2levels"))]
    names(arglist)[names(arglist) == "blocks_attr"] <- "attr"
    arguments$constraints$blocks <- InitErgmConstraint.blocks(nw, arglist)
  }
  if(is.null(arguments$constraints$blocks)) {
    arguments$constraints$blocks <- InitErgmConstraint.blocks(nw, list(attr = trim_env(~0)))
  }

  # check for old name
  if(is.null(arguments$constraints$strat) && !is.null(arguments$constraints$Strat)) {
    arguments$constraints$strat <- arguments$constraints$Strat
  }

  # if strat has not already been initialized, or if related arguments are passed directly to the proposal, (re)initialize it now
  if(length(intersect(names(arguments), c("strat_attr", "pmat", "empirical"))) > 0) {
    arglist <- arguments[intersect(names(arguments), c("strat_attr", "pmat", "empirical"))]
    names(arglist)[names(arglist) == "strat_attr"] <- "attr"
    arguments$constraints$strat <- InitErgmConstraint.strat(nw, arglist)
  }
  if(is.null(arguments$constraints$strat)) {
    arguments$constraints$strat <- InitErgmConstraint.strat(nw, list(attr = trim_env(~0)))
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
                   strattailattrs = as.integer(arguments$constraints$strat$tailattrs - 1L),
                   stratheadattrs = as.integer(arguments$constraints$strat$headattrs - 1L),
                   probvec = as.double(arguments$constraints$strat$probvec),
                   nattrcodes = as.integer(arguments$constraints$strat$nlevels),
                   strat_vattr = as.integer(arguments$constraints$strat$nodecov - 1L),
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
                   amat = as.integer(t(pairs_mat)),
                   dyad_indep = as.integer(dyad_indep))

  proposal
}

#' @templateVar name CondDegree
#' @aliases InitErgmProposal.CondDegree
#' @title MHp for degree constraints
#' @description MHp for \eqn{constraints= ~degree}. Propose either 4 toggles (MH_CondDegreeTetrad) or 6 toggles
#'   (MH_CondDegreeHexad) at once. For undirected networks, propose 4 toggles (MH_CondDegreeTetrad).
#'   MH_CondDegreeTetrad selects two edges with no nodes in common, A1-A2 and B1-B2, s.t. A1-B2 and B1-A2 are
#'   not edges, and propose to replace the former two by the latter two. MH_CondDegreeHexad selects three edges
#'   A1->A2, B1->B2, C1->C2 at random and rotate them to A1->B2, B1->C2, and C1->A2.
#' @template ergmProposal-general
NULL
InitErgmProposal.CondDegree <- function(arguments, nw) {
  proposal <- list(name = "CondDegree", inputs=NULL)
  proposal
}

#' @templateVar name CondDegreeMix
#' @aliases InitErgmProposal.CondDegreeMix
#' @title MHp for degree mix constraints
#' @description MHp for \eqn{constraints= ~degreesmix}. Similar to `InitErgmProposal.CondDegree`, except that
#'   the toggle is proposed only if the mixing matrix of degrees is preserved before and after the toggle.
#' @template ergmProposal-general
NULL
InitErgmProposal.CondDegreeMix <- function(arguments, nw) {
  proposal <- list(name = "CondDegreeMix",
    inputs=get.vertex.attribute(nw,arguments$constraints$degreesmix$attrib))
  proposal
}

#' @templateVar name CondOutDegree
#' @aliases InitErgmProposal.CondOutDegree
#' @title MHp for odegree constraints
#' @description MHp for \eqn{constraints= ~odegrees}. For directed networks, randomly select two dyads with a
#'   common tail node, one having an edge and one not, and propose to swap the tie from one head to the other.
#' @template ergmProposal-general
NULL
InitErgmProposal.CondOutDegree <- function(arguments, nw) {
  proposal <- list(name = "CondOutDegree", inputs=NULL)
  if (!is.directed(nw)) # Really, this should never trigger, since the InitErgmConstraint function should check.
    ergm_Init_abort("The CondOutDegree proposal function does not work with an",
          "undirected network.")
  proposal
}

#' @templateVar name CondInDegree
#' @aliases InitErgmProposal.CondInDegree
#' @title MHp for idegree constraints
#' @description MHp for \eqn{constraints= ~idegrees}. For directed networks, randomly select two dyads with a
#'   common head node, one having an edge one not, and propose to swap the tie from one tail to the other.
#' @template ergmProposal-general
NULL
InitErgmProposal.CondInDegree <- function(arguments, nw) {
  proposal <- list(name = "CondInDegree", inputs=NULL)
  if (!is.directed(nw)) # Really, this should never trigger, since the InitErgmConstraint function should check.
    ergm_Init_abort("The CondInDegree proposal function does not work with an",
          "undirected network.")
  proposal
}

#' @templateVar name CondB1Degree
#' @aliases InitErgmProposal.CondB1Degree
#' @title MHp for b1degree constraints
#' @description MHp for \eqn{constraints= ~b1degrees}. For bipartite networks, randomly select an edge {B1i, B2j}
#'   and an empty dyad with the same node B1i, {B1i, B2k}, and propose to toggle both {B1i, B2j} and {B1i, B2k}.
#'   This ensures that the degrees of individual nodes in mode 1 are preserved.
#' @template ergmProposal-general
NULL
InitErgmProposal.CondB1Degree <- function(arguments, nw) {
  proposal <- list(name = "CondB1Degree", inputs=NULL)
  if (!is.bipartite(nw)) # Really, this should never trigger, since the InitErgmConstraint function should check.
    ergm_Init_abort("The CondB1Degree proposal function does not work with a non-bipartite network.")
  
  proposal
}

#' @templateVar name CondB2Degree
#' @aliases InitErgmProposal.CondB2Degree
#' @title MHp for b2degree constraints
#' @description MHp for \eqn{constraints= ~b2degrees}. For bipartite networks, randomly select an edge {B1j, B2i}
#'   and an empty dyad with the same node B2i, {B1k, B2i}, and propose to toggle both {B1j, B2i} and {B1k, B2i}.
#'   This ensures that the degrees of individual nodes in mode 2 are preserved.
#' @template ergmProposal-general
NULL
InitErgmProposal.CondB2Degree <- function(arguments, nw) {
  proposal <- list(name = "CondB2Degree", inputs=NULL)
  if (!is.bipartite(nw)) # Really, this should never trigger, since the InitErgmConstraint function should check.
    ergm_Init_abort("The CondB2Degree proposal function does not work with a non-bipartite network.")
  proposal
}

#' @templateVar name CondDegreeDist
#' @aliases InitErgmProposal.CondDegreeDist
#' @title MHp for degreedist constraints
#' @description MHp for \eqn{constraints= ~degreedist}. Randomly select a node (T) and its edge (E).  If the head
#'   node of the edge (H) has 1 degree more than another randomly select node (A), and A is disconnected to both
#'   T and H, then propose to toggle E and the dyad between T and A.
#' @template ergmProposal-general
NULL
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

#' @templateVar name CondInDegreeDist
#' @aliases InitErgmProposal.CondInDegreeDist
#' @title MHp for idegreedist constraints
#' @description MHp for \eqn{constraints= ~idegreedist}. For directed networks, similar to
#'   `InitErgmProposal.CondDegreeDist`, except for indegree case
#' @template ergmProposal-general
NULL
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

#' @templateVar name CondOutDegreeDist
#' @aliases InitErgmProposal.CondOutDegreeDist
#' @title MHp for odegreedist constraints
#' @description MHp for \eqn{constraints= ~odegreedist}. For directed networks, similar to
#'   `InitErgmProposal.CondDegreeDist`, except for outdegree case
#' @template ergmProposal-general
NULL
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

#' @templateVar name ConstantEdges
#' @aliases InitErgmProposal.ConstantEdges
#' @title MHp for edges constraints
#' @description MHp for \eqn{constraints= ~edges}. Propose pairs of toggles that keep number of edges the same.
#'   This is done by:
#'   a. choosing an existing edge at random;
#'   b. repeatedly choosing dyads at random until one is found that does not have an edge; and
#'   c. proposing toggling both these dyads. Note that step b. will be very inefficient if the network is nearly
#'      complete, so this proposal is NOT recommended for such networks. However, most network datasets are
#'      sparse, so this is not likely to be an issue.
#' @template ergmProposal-general
NULL
InitErgmProposal.ConstantEdges <- function(arguments, nw) {
  proposal <- list(name = "ConstantEdges", bd = ergm_bd_init(arguments, nw))
  proposal
}

#' @templateVar name HammingConstantEdges
#' @aliases InitErgmProposal.HammingConstantEdges
#' @title TODO
#' @description TODO
#' @template ergmProposal-general
NULL
InitErgmProposal.HammingConstantEdges <- function(arguments, nw) {
  proposal <- list(name = "HammingConstantEdges", inputs=NULL)
  if(is.bipartite(nw)){
    proposal$name <- "BipartiteHammingConstantEdges"
  }
  proposal
}

#' @templateVar name HammingTNT
#' @aliases InitErgmProposal.HammingTNT
#' @title TODO
#' @description TODO
#' @template ergmProposal-general
NULL
InitErgmProposal.HammingTNT <- function(arguments, nw) {
  proposal <- list(name = "HammingTNT", inputs=NULL)
  if(is.bipartite(nw)){
    proposal$name <- "BipartiteHammingTNT"
  }
  proposal
}
