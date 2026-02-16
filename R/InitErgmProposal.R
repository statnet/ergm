#  File R/InitErgmProposal.R in package ergm, part of the Statnet suite of
#  packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
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

#' @templateVar name randomtoggle
#' @aliases InitErgmProposal.randomtoggle
#' @title Propose a randomly selected dyad to toggle
#' @description Propose a randomly selected dyad to toggle
#' @template ergmProposal-general
NULL
InitErgmProposal.randomtoggle <- function(arguments, nw){
  c(list(name = "randomtoggle", dyadgen = ergm_dyadgen_select(arguments, nw), bd = ergm_bd_init(arguments, nw)),
    ergm_constrain_changestats(arguments))
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
  c(list(name = "TNT", dyadgen = ergm_dyadgen_select(arguments, nw), bd = ergm_bd_init(arguments, nw)),
    ergm_constrain_changestats(arguments))
}

#' @templateVar name BDStratTNT
#' @aliases InitErgmProposal.BDStratTNT
#' @title TNT proposal with degree bounds, stratification, and a blocks constraint
#' @description Implements a TNT proposal with any subset of the following features:
#'   1. upper bounds on degree, specified via the [`bd`][bd-ergmConstraint]
#'      constraint's `maxout`, `maxin`, and `attribs` arguments;
#'   2. stratification of proposals according to mixing type on a vertex attribute,
#'      specified via the [`strat`][strat-ergmHint] hint;
#'   3. fixation of specified mixing types on a(nother) vertex attribute, specified
#'      via the [`blocks`][blocks-ergmConstraint] constraint.
#' @template ergmProposal-general
NULL
InitErgmProposal.BDStratTNT <- function(arguments, nw) {
  ## constraints are dyad-independent if bd is not being used
  dyad_indep <- is.null(arguments$constraints$bd)

  ## handle defaults for hints
  NVL(arguments$constraints$bd) <- InitErgmConstraint.bd(nw, list())
  NVL(arguments$constraints$blocks) <- InitErgmConstraint.blocks(nw, list(attr = trim_env(~0)))
  NVL(arguments$constraints$strat) <- InitErgmConstraint.strat(nw, list(attr = trim_env(~0)))

  attribs <- NVL(arguments$constraints$bd$attribs,
                 matrix(TRUE, ncol = 1L, nrow = network.size(nw)))

  maxout <- NVL(arguments$constraints$bd$maxout, network.size(nw))
  maxout %[f]% is.na <- network.size(nw)
  maxout <- matrix(rep(maxout, length.out = length(attribs)), ncol = ncol(attribs))

  maxin <- NVL(arguments$constraints$bd$maxin, network.size(nw))
  maxin %[f]% is.na <- network.size(nw)
  maxin <- matrix(rep(maxin, length.out = length(attribs)), ncol = ncol(attribs))

  bd_vattr <- which(attribs, arr.ind = TRUE)
  bd_vattr <- bd_vattr[order(bd_vattr[, 1L]), 2L]
  bd_nlevels <- NCOL(attribs)

  ## allowed bd pairings
  bd_mixmat <- matrix(TRUE, nrow = NCOL(attribs), ncol = NCOL(attribs))
  if(!is.directed(nw)) {
    bd_mixmat[lower.tri(bd_mixmat, diag = FALSE)] <- FALSE
  }

  bd_pairs <- which(bd_mixmat, arr.ind = TRUE)
  bd_tails <- bd_pairs[, 1L]
  bd_heads <- bd_pairs[, 2L]

  ## need to handle undirected case as for blocks, but bipartite along with unip for now
  bd_offdiag_pairs <- which(bd_tails != bd_heads)

  bd_allowed_tails <- c(bd_tails, if(!is.directed(nw)) bd_heads[bd_offdiag_pairs])
  bd_allowed_heads <- c(bd_heads, if(!is.directed(nw)) bd_tails[bd_offdiag_pairs])

  ## number of bd mixtypes that need to be considered when strat and blocks
  ## mixing types are off-diag and on-diag, respectively
  bd_nmixtypes <- c(length(bd_allowed_tails), length(bd_tails))

  blocks_vattr <- arguments$constraints$blocks$nodecov
  amat <- arguments$constraints$blocks$amat

  blocks_pairs <- which(amat, arr.ind = TRUE)
  if(!is.directed(nw) && !is.bipartite(nw)) {
    blocks_pairs <- blocks_pairs[blocks_pairs[, 1L] <= blocks_pairs[, 2L], , drop = FALSE]
  }
  blocks_tails <- blocks_pairs[, 1L]
  blocks_heads <- blocks_pairs[, 2L]

  blocks_nlevels <- NROW(amat)
  blocks_node_counts <- tabulate(blocks_vattr, nbins = blocks_nlevels)

  pairs_to_keep <- (blocks_tails != blocks_heads
                    & blocks_node_counts[blocks_tails] > 0
                    & blocks_node_counts[blocks_heads] > 0) |
                   (blocks_tails == blocks_heads
                    & blocks_node_counts[blocks_tails] > 1)
  blocks_tails <- blocks_tails[pairs_to_keep]
  blocks_heads <- blocks_heads[pairs_to_keep]

  blocks_offdiag_pairs <- which(blocks_tails != blocks_heads)

  blocks_allowed_tails <- c(blocks_tails,
                            if(!is.bipartite(nw) && !is.directed(nw)) blocks_heads[blocks_offdiag_pairs])
  blocks_allowed_heads <- c(blocks_heads,
                            if(!is.bipartite(nw) && !is.directed(nw)) blocks_tails[blocks_offdiag_pairs])

  ## number of blocks mixtypes that need to be considered
  ## when strat mixing type is off-diag and on-diag, respectively
  blocks_nmixtypes <- c(length(blocks_allowed_tails), length(blocks_tails))

  strat_vattr <- arguments$constraints$strat$nodecov
  strat_nlevels <- arguments$constraints$strat$nlevels

  strat_vattr <- strat_vattr - 1L
  blocks_vattr <- blocks_vattr - 1L
  bd_vattr <- bd_vattr - 1L

  combined_vattr <- strat_vattr*blocks_nlevels*bd_nlevels + blocks_vattr*bd_nlevels + bd_vattr
  combined_nlevels <- strat_nlevels*blocks_nlevels*bd_nlevels
  combined_vattr_counts <- tabulate(combined_vattr + 1L, nbins = combined_nlevels)  

  proposal <- list(name = "BDStratTNT",
                   inputs = NULL, # passed by name below
                   strat_vattr = as.integer(strat_vattr),
                   strat_nlevels = as.integer(strat_nlevels),
                   strat_tails = as.integer(arguments$constraints$strat$tailattrs - 1L),
                   strat_heads = as.integer(arguments$constraints$strat$headattrs - 1L),
                   blocks_vattr = as.integer(blocks_vattr),
                   blocks_nlevels = as.integer(blocks_nlevels),
                   blocks_tails = as.integer(blocks_allowed_tails - 1L),
                   blocks_heads = as.integer(blocks_allowed_heads - 1L),
                   blocks_nmixtypes = as.integer(blocks_nmixtypes),
                   bd_vattr = as.integer(bd_vattr),
                   bd_nlevels = as.integer(bd_nlevels),
                   bd_tails = as.integer(bd_allowed_tails - 1L),
                   bd_heads = as.integer(bd_allowed_heads - 1L),
                   bd_nmixtypes = as.integer(bd_nmixtypes),
                   combined_vattr = as.integer(combined_vattr),
                   combined_nlevels = as.integer(combined_nlevels),
                   combined_vattr_counts = as.integer(combined_vattr_counts),
                   maxout = as.integer(maxout),
                   maxin = as.integer(maxin),
                   probvec = as.double(arguments$constraints$strat$probvec),
                   empirical_flag = as.integer(arguments$constraints$strat$empirical),
                   amat = as.integer(t(amat)),
                   dyad_indep = as.integer(dyad_indep))

  proposal
}

#' @templateVar name CondDegree
#' @aliases InitErgmProposal.CondDegree
#' @title MHp for degree constraints
#' @description MHp for `constraints= ~degree`. Propose either 4 toggles (MH_CondDegreeTetrad) or 6 toggles
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
#' @description MHp for `constraints= ~degreesmix`. Similar to `InitErgmProposal.CondDegree`, except that
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
#' @description MHp for `constraints= ~odegrees`. For directed networks, randomly select two dyads with a
#'   common tail node, one having an edge and one not, and propose to swap the tie from one head to the other.
#' @template ergmProposal-general
NULL
InitErgmProposal.CondOutDegree <- function(arguments, nw) {
  proposal <- list(name = "CondOutDegree", inputs=NULL)
  if (!is.directed(nw)) # Really, this should never trigger, since the InitErgmConstraint function should check.
    ergm_Init_stop("The CondOutDegree proposal function does not work with an",
          "undirected network.")
  proposal
}

#' @templateVar name CondInDegree
#' @aliases InitErgmProposal.CondInDegree
#' @title MHp for idegree constraints
#' @description MHp for `constraints= ~idegrees`. For directed networks, randomly select two dyads with a
#'   common head node, one having an edge one not, and propose to swap the tie from one tail to the other.
#' @template ergmProposal-general
NULL
InitErgmProposal.CondInDegree <- function(arguments, nw) {
  proposal <- list(name = "CondInDegree", inputs=NULL)
  if (!is.directed(nw)) # Really, this should never trigger, since the InitErgmConstraint function should check.
    ergm_Init_stop("The CondInDegree proposal function does not work with an",
          "undirected network.")
  proposal
}

#' @templateVar name CondB1Degree
#' @aliases InitErgmProposal.CondB1Degree
#' @title MHp for b1degree constraints
#' @description MHp for `constraints= ~b1degrees`. For bipartite networks, randomly select an edge \eqn{(B_{1i},B_{2j})}
#'   and an empty dyad with the same node B1i, \eqn{(B_{1i},B_{2k})}, and propose to toggle both \eqn{(B_{1i},B_{2j})} and \eqn{(B_{1i},B_{2k})}.
#'   This ensures that the degrees of individual nodes in mode 1 are preserved.
#' @template ergmProposal-general
NULL
InitErgmProposal.CondB1Degree <- function(arguments, nw) {
  proposal <- list(name = "CondB1Degree", inputs=NULL)
  if (!is.bipartite(nw)) # Really, this should never trigger, since the InitErgmConstraint function should check.
    ergm_Init_stop("The CondB1Degree proposal function does not work with a non-bipartite network.")
  
  proposal
}

#' @templateVar name CondB2Degree
#' @aliases InitErgmProposal.CondB2Degree
#' @title MHp for b2degree constraints
#' @description MHp for `constraints= ~b2degrees`. For bipartite networks, randomly select an edge \eqn{(B_{1j},B_{2i})}
#'   and an empty dyad with the same node B2i, \eqn{(B_{1k},B_{2i})}, and propose to toggle both \eqn{(B_{1j},B_{2i})} and \eqn{(B_{1k},B_{2i})}.
#'   This ensures that the degrees of individual nodes in mode 2 are preserved.
#' @template ergmProposal-general
NULL
InitErgmProposal.CondB2Degree <- function(arguments, nw) {
  proposal <- list(name = "CondB2Degree", inputs=NULL)
  if (!is.bipartite(nw)) # Really, this should never trigger, since the InitErgmConstraint function should check.
    ergm_Init_stop("The CondB2Degree proposal function does not work with a non-bipartite network.")
  proposal
}

#' @templateVar name CondDegreeDist
#' @aliases InitErgmProposal.CondDegreeDist
#' @title MHp for degreedist constraints
#' @description MHp for `constraints= ~degreedist`. Randomly select a node (T) and its edge (E).  If the head
#'   node of the edge (H) has 1 degree more than another randomly select node (A), and A is disconnected to both
#'   T and H, then propose to toggle E and the dyad between T and A.
#' @template ergmProposal-general
NULL
InitErgmProposal.CondDegreeDist <- function(arguments, nw) {
  proposal <- list(name = "CondDegreeDist", inputs=NULL)
  if (is.directed(nw)) {
    ergm_Init_warning("Using the 'degreedist' constraint with a directed network ",
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
#' @description MHp for `constraints= ~idegreedist`. For directed networks, similar to
#'   `InitErgmProposal.CondDegreeDist`, except for indegree case
#' @template ergmProposal-general
NULL
InitErgmProposal.CondInDegreeDist <- function(arguments, nw) {
  proposal <- list(name = "CondInDegreeDist", inputs=NULL)
  if (!is.directed(nw)) {
    ergm_Init_warning("Using the 'idegreedist' constraint with an undirected network ",
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
#' @description MHp for `constraints= ~odegreedist`. For directed networks, similar to
#'   `InitErgmProposal.CondDegreeDist`, except for outdegree case
#' @template ergmProposal-general
NULL
InitErgmProposal.CondOutDegreeDist <- function(arguments, nw) {
  proposal <- list(name = "CondOutDegreeDist", inputs=NULL)
  if (!is.directed(nw)) {
    ergm_Init_warning("Using the 'odegreedist' constraint with an undirected network n",
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
#' @description MHp for `constraints= ~edges`. Propose pairs of toggles that keep number of edges the same.
#'   This is done by:
#'
#'   1. choosing an existing edge at random;
#'   2. choosing a dyad at random that does not have an edge; and
#'   3. proposing toggling both these dyads.
#' @template ergmProposal-general
NULL
InitErgmProposal.ConstantEdges <- function(nw, arguments, ...) {
  list(name = "ConstantEdges", dyadgen = ergm_dyadgen_select(arguments, nw), bd = ergm_bd_init(arguments, nw))
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

#' @templateVar name SPDyad
#' @aliases InitErgmProposal.SPDyad
#' @title A proposal alternating between TNT and a triad-focused
#'   proposal
#' @description The specified proportion of the time, the proposal
#'   proceeds along the lines of \insertCite{WaAt13a;textual}{ergm},
#'   albeit with different weighting. A dyad is selected uniformly at
#'   random from among those dyads with at least one shared
#'   partnership or transitivity of the specified type. This is likely
#'   to be more efficient for a model with excess triangles.
#'
#' @references \insertAllCited{}
#'
#' @template ergmProposal-general
NULL
InitErgmProposal.SPDyad <- function(arguments, nw) {
  c(list(name = "SPDyad",
         inputs = c(NVL(arguments$constraints$triadic$triFocus,0.25)),
         iinputs = SPTYPE_CODE[arguments$constraints$triadic$type],
         dyadgen = ergm_dyadgen_select(arguments, nw), bd = ergm_bd_init(arguments, nw)),
    ergm_constrain_changestats(arguments, .spcache.aux(arguments$constraints$triadic$type)))
}
