#  File R/obs.constraints.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################

.delete_from_conform_rhs <- function(f, del){
  ff <- filter_rhs.formula(f, `!=`, del)
  if(length(ff)!=length(f)) ff[[length(ff)+1]] <- as.name('.')
  ff
}

.handle.auto.constraints <- function(nw,
                                     constraints=~.,
                                     obs.constraints=~.-observed,
                                     target.stats=NULL) {

  # We have constraint information.
  constraints <- NVL3(nw%ergmlhs%"constraints", nonsimp_update.formula(., constraints), constraints)

  if(!is.null(obs.constraints)){
    # We have observational process information.
    obs.constraints <- NVL3(nw%ergmlhs%"obs.constraints", nonsimp_update.formula(., obs.constraints), obs.constraints)
    
    # Observation process handling only needs to happen if the
    # sufficient statistics are not specified. If the sufficient
    # statistics are specified, the nw's dyad states are irrelevant.
    if(!is.null(target.stats)){
      if(obs.constraints!=~. || network.naedgecount(nw)) message("Target statistics specified in a network with missing dyads and/or a nontrivial observation process. Since (by sufficiency) target statistics provide all the information needed to fit the model, missingness and observation process will not affect estimation.")
      if(network.naedgecount(nw)) nw[as.matrix(is.na(nw),matrix.type="edgelist")] <- 0
      obs.constraints <- ~.
    }

    # If no missing edges, remove the "observed" constraint.
    if(network.naedgecount(nw)==0){
      obs.constraints <- .delete_from_conform_rhs(obs.constraints, "observed")
    }
    
    constraints.obs<-obs.constraints
    ult(constraints.obs) <- call('+', ult(constraints.obs), ult(constraints))
    constraints.obs <- .delete_from_conform_rhs(constraints.obs, ".")
    if(constraints==constraints.obs) constraints.obs<-NULL
    
  }else constraints.obs<-NULL
  
  list(nw = nw, constraints = constraints, constraints.obs = constraints.obs)
}

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
