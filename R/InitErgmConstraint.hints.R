#  File R/InitErgmConstraint.hints.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################
InitErgmConstraint.TNT<-function(nw, arglist, ...){
  .Deprecate_once("sparse")
  InitErgmConstraint.sparse(nw, arglist, ...)
}

InitErgmConstraint.sparse<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist)
  list(dependence = FALSE, priority=10, impliedby=c("sparse", "edges", "degrees", "edges", "idegrees", "odegrees", "b1degrees", "b2degrees", "idegreedist", "odegreedist", "degreedist", "b1degreedist", "b2degreedist"), constrain="sparse")
}

InitErgmConstraint.Strat<-function(nw, arglist, ...){
  .Deprecate_once("strat")
  InitErgmConstraint.strat(nw, arglist, ...)
}

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
    ergm_Init_abort(sQuote("pmat"), " does not have the correct dimensions for ", sQuote("Strat_attr"), ".")    
  }
  
  if(!is.bipartite(nw) && !is.directed(nw)) {
    # for undirected unipartite, symmetrize pmat and then set the sub-diagonal to zero
    pmat <- (pmat + t(pmat))/2
    pmat[lower.tri(pmat)] <- 0
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
  
  list(dependence = FALSE, 
       priority = 10, 
       tailattrs = tailattrs,
       headattrs = headattrs,
       probvec = probvec,
       nlevels = nlevels,
       nodecov = strat_nodecov,
       empirical = empirical,
       constrain = "strat")
}
