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
    strat_row_nodecov <- NVL2(attr, ergm_get_vattr(attr, nw, bip="b1"), rep(1, nw %n% "bipartite"))
    strat_col_nodecov <- NVL2(attr, ergm_get_vattr(attr, nw, bip="b2"), rep(1, network.size(nw) - (nw %n% "bipartite")))

    strat_row_levels <- sort(unique(strat_row_nodecov))
    strat_col_levels <- sort(unique(strat_col_nodecov))

    strat_levels <- c(strat_row_levels, strat_col_levels)

    strat_row_nodecov <- match(strat_row_nodecov, strat_row_levels)
    strat_col_nodecov <- match(strat_col_nodecov, strat_col_levels)
    strat_nodecov <- c(strat_row_nodecov, strat_col_nodecov + length(strat_row_levels))
  } else {
    strat_row_nodecov <- NVL2(attr, ergm_get_vattr(attr, nw), rep(1, network.size(nw)))
    strat_col_nodecov <- strat_row_nodecov

    strat_row_levels <- sort(unique(strat_row_nodecov))
    strat_col_levels <- strat_row_levels

    strat_levels <- strat_row_levels
    
    strat_row_nodecov <- match(strat_row_nodecov, strat_row_levels)
    strat_col_nodecov <- strat_row_nodecov
    strat_nodecov <- strat_row_nodecov
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

  list(dependence = FALSE, 
       priority = 10, 
       nmixtypes = nmixingtypes,
       tailattrs = tailattrs,
       headattrs = headattrs,
       probvec = probvec,
       nlevels = ncodes,
       nodecountsbycode = nodecountsbycode,
       nodeindicesbycode = nodeindicesbycode,
       nodecov = strat_nodecov,
       indmat = indmat,
       empirical = empirical,
       constrain = "strat")
}
