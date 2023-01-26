#  File R/obs.constraints.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2023 Statnet Commons
################################################################################

.handle.auto.constraints <- function(nw,
                                     constraints=trim_env(~.),
                                     obs.constraints=trim_env(~.-observed),
                                     target.stats=NULL,
                                     default.dot=c("first", "last", "none")){
  default.dot <- match.arg(default.dot)

  .preproc_constraints <- function(...){
    # Embed the LHS as a .select() constraint.
    tll <- list(...) %>% compact() %>% map(.embed_constraint_lhs) %>% map(list_rhs.formula)
    # If no missing edges, remove the "observed" constraint.
    if(network.naedgecount(nw)==0) tll <- map(tll, .delete_term, "observed")

    if(length(tll) == 0) return(NULL)

    # Go through the constraint lists, substituting the earlier ones into the dots in the later ones.
    otl <- tll[[1]]
    for(i in seq_along(tll)[-1]){
      ntl <- tll[[i]]

      # Find the substitution positions.
      pos <- which(ntl %>% map_chr(~as.character(.)[1]) == ".")

      # Don't substitute at negative dots.
      pos_sign <- attr(ntl, "sign")[pos]
      pos <- pos[pos_sign>0]

      # If no dots, use default behaviour unless there is a negative dot.
      if(!length(pos) && all(pos_sign>0))
        ntl <- switch(default.dot,
                     first = c(otl, ntl),
                     last = c(ntl, otl),
                     none = ntl)

      # Substitute (working backwards, to prevent pos from changing).
      for(p in rev(pos))
        ntl <- c(ntl[seq_len(pos-1)], otl, ntl[-seq_len(pos)])

      otl <- ntl
    }

    # Delete remaining dots (including the negative ones).
    .delete_term(otl, ".")
  }

  tl <- .preproc_constraints(nw%ergmlhs%"constraints", constraints)
  obs.tl <- .preproc_constraints(nw%ergmlhs%"obs.constraints", obs.constraints)

  # Do any of the observational constraints formulas have terms?
  if(length(obs.tl)){
    # Observation process handling only needs to happen if the
    # sufficient statistics are not specified. If the sufficient
    # statistics are specified, the nw's dyad states are irrelevant.
    if(!is.null(target.stats)){
      message("Target statistics specified in a network with missing dyads and/or a nontrivial observation process. Since (by sufficiency) target statistics provide all the information needed to fit the model, missingness and observation process will not affect estimation.")
      if(network.naedgecount(nw)) nw[as.matrix(is.na(nw),matrix.type="edgelist")] <- 0
      obs.tl <- NULL
    }
  }

  list(nw = nw, conterms = tl, conterms.obs = if(length(obs.tl)) c(obs.tl, tl))
}

has.obs.constraints <- function(...) length(.handle.auto.constraints(...)$conterms.obs) > 0

.align.target.stats.offset <- function(model, target.stats){
  om <- model$etamap$offsetmap
  cno <- param_names(model, canonical=TRUE)
  cn <- param_names(model, canonical=TRUE)[!om]
  target.stats <- na.omit(vector.namesmatch(target.stats, cn))
  tmp <- rep(NA, length(om))
  tmp[!om] <- target.stats
  names(tmp) <- cno
  tmp
}
