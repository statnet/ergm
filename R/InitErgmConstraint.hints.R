#  File R/InitErgmConstraint.hints.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################
InitErgmConstraint.TNT<-function(nw, arglist, ...){
  .Deprecate_once("sparse")
  InitErgmConstraint.sparse(nw, arglist, ...)
}

#' @templateVar name sparse
#' @title Sparse network
#' @description The network is sparse. This typically results in a Tie-Non-Tie (TNT) proposal regime.
#'
#' @usage
#' # sparse
#'
#' @template ergmHint-general
#'
#' @concept dyad-independent
InitErgmConstraint.sparse<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist)
  list(priority=10, impliedby=c("sparse", "edges", "degrees", "edges", "idegrees", "odegrees", "b1degrees", "b2degrees", "idegreedist", "odegreedist", "degreedist", "b1degreedist", "b2degreedist"), constrain="sparse")
}

#' @templateVar name triadic
#' @title Network with strong clustering (triad-closure) effects
#' @description The network has a high clustering coefficient. This typically results in alternating between the Tie-Non-Tie (TNT) proposal and a triad-focused proposal along the lines of that of \insertCite{WaAt13a;textual}{ergm}.
#'
#' @usage
#' # triadic(triFocus = 0.25, type="OTP")
#'
#' @param triFocus A number between 0 and 1, indicating how often triad-focused proposals should be made relative to the standard proposals.
#' @template ergmTerm-sp-type
#'
#' @template ergmTerm-sp-types
#'
#' @template ergmHint-general
#'
#' @references \insertAllCited{}
#'
#' @concept dyad-dependent
InitErgmConstraint.triadic<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist, bipartite = FALSE,
                      varnames = c("triFocus", "type"),
                      vartypes = c("numeric", "character"),
                      defaultvalues = list(0.25, "OTP"),
                      required = c(FALSE, FALSE))
  if(!is.directed(nw)) a$type <- "UTP"
  list(triFocus=a$triFocus, type = a$type, priority=4, constrain="triadic")
}

#' @templateVar name triadic
#' @template ergmHint-rdname
#' @aliases .triadic-ergmHint
#' @usage
#' # .triadic(triFocus = 0.25, type = "OTP")
#' @section `.triadic()` versus `triadic()`: If given a bipartite
#'   network, the dotted form will skip silently, whereas the plain
#'   form will raise an error, since triadic effects are not possible
#'   in bipartite networks. The dotted form is thus suitable as a
#'   default argument when the bipartitedness of the network is not
#'   known *a priori*.
InitErgmConstraint..triadic<-function(nw, arglist, ...){
  if(is.bipartite(nw)) NULL
  else InitErgmConstraint.triadic(nw, arglist, ...)
}

InitErgmConstraint.Strat<-function(nw, arglist, ...){
  .Deprecate_once("strat")
  InitErgmConstraint.strat(nw, arglist, ...)
}

#' @templateVar name strat
#' @title Stratify Proposed Toggles by Mixing Type on a Vertex Attribute
#' @description Proposed toggles are stratified according to mixing type
#'              on a vertex attribute.
#'
#' @details The user may pass a vertex attribute `attr` as an argument
#'   (the default for `attr` gives every vertex the same attribute
#'   value), and may also pass a matrix of weights `pmat` (the default
#'   for `pmat` gives equal weight to each mixing type). See
#'   [Specifying Vertex Attributes and Levels for
#'   details][nodal_attributes] on specifying vertex attributes. The
#'   matrix `pmat`, if specified, must have the same dimensions as a
#'   mixing matrix for the network and vertex attribute under 
#'   consideration, and the correspondence between rows and columns of
#'   `pmat` and values of `attr` is the same as for a mixing matrix.
#'
#'   The interpretation is that `pmat[i,j]/sum(pmat)` is the probability of
#'   proposing a toggle for mixing type `(i,j)`. (For undirected, unipartite
#'   networks, `pmat` is first symmetrized, and then entries below the diagonal
#'   are set to zero. Only entries on or above the diagonal of the symmetrized
#'   `pmat` are considered when making proposals. This accounts for the
#'   convention that mixing is undirected in an undirected, unipartite network:
#'   a tail of type `i` and a head of type `j` has the same mixing type
#'   as a tail of type `j` and a head of type `i`.)
#'
#'   As an alternative way of specifying `pmat`, the user may pass
#'   `empirical = TRUE` to use the mixing matrix of the network beginning
#'   the MCMC chain as `pmat`. In order for this to work, that network should
#'   have a reasonable (in particular, nonempty) edge set.
#'
#'   While some mixing types may be assigned zero proposal probability
#'   (either with a direct specification of `pmat` or with `empirical = TRUE`),
#'   this will not be recognized as a constraint by all components of `ergm`,
#'   and should be used with caution.
#'
#' @usage
#' # strat(attr=NULL, pmat=NULL, empirical=FALSE)
#'
#' @template ergmHint-general
#'
#' @concept dyad-independent
InitErgmConstraint.strat <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attr", "pmat", "empirical"),
                      vartypes = c(ERGM_VATTR_SPEC, "matrix,table", "logical"),
                      defaultvalues = list(NULL, NULL, FALSE),
                      required = c(TRUE, FALSE, FALSE))
  attr <- a$attr; pmat <- a$pmat; empirical <- a$empirical

  if(is.bipartite(nw)) {
    strat_row_nodecov <- ergm_get_vattr(attr, nw, bip = "b1")
    strat_col_nodecov <- ergm_get_vattr(attr, nw, bip = "b2")

    strat_row_levels <- sort(unique(strat_row_nodecov))
    strat_col_levels <- sort(unique(strat_col_nodecov))

    strat_nodecov <- c(match(strat_row_nodecov, strat_row_levels),
                       match(strat_col_nodecov, strat_col_levels) + length(strat_row_levels))

    strat_levels <- c(strat_row_levels, strat_col_levels)
  } else {
    strat_nodecov <- ergm_get_vattr(attr, nw)

    strat_levels <- sort(unique(strat_nodecov))

    strat_nodecov <- match(strat_nodecov, strat_levels)

    strat_row_levels <- strat_levels
    strat_col_levels <- strat_levels
  }

  pmat <- NVL(pmat, matrix(1, nrow = length(strat_row_levels), ncol = length(strat_col_levels)))

  if(NROW(pmat) != length(strat_row_levels) || NCOL(pmat) != length(strat_col_levels)) {
    ergm_Init_stop(sQuote("pmat"), " does not have the correct dimensions for ", sQuote("attr"), ".")
  }

  if(!is.bipartite(nw) && !is.directed(nw)) {
    # for undirected unipartite, symmetrize pmat and then set the sub-diagonal to zero
    pmat <- (pmat + t(pmat))/2
    pmat[lower.tri(pmat, diag = FALSE)] <- 0
  }

  # record the tail and head attr code for each mixing type with positive probability
  prob_inds <- which(pmat > 0, arr.ind = TRUE)
  tailattrs <- prob_inds[,1]
  headattrs <- prob_inds[,2]
  probvec <- pmat[prob_inds]

  if(is.bipartite(nw)) {
    headattrs <- headattrs + length(strat_row_levels)
  }

  # record the number of unique attr codes
  nlevels <- length(strat_levels)

  list(priority = 4,
       tailattrs = tailattrs,
       headattrs = headattrs,
       probvec = probvec,
       nlevels = nlevels,
       nodecov = strat_nodecov,
       empirical = empirical,
       constrain = "strat")
}
