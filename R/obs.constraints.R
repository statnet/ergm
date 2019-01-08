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
  if("constraints" %in% list.network.attributes(nw)){
    constraints <- nonsimp_update.formula(nw %n% "constraints", constraints)
  }

  if(!is.null(obs.constraints)){
    # We have observational process information.
    if("obs.constraints" %in% list.network.attributes(nw)){
      obs.constraints <- nonsimp_update.formula(nw %n% "obs.constraints", obs.constraints)
    }
    
    # Observation process handling only needs to happen if the
    # sufficient statistics are not specified. If the sufficient
    # statistics are specified, the nw's dyad states are irrelevant.
    if(!is.null(target.stats)){
      if(network.naedgecount(nw)){
        warning("Target statistics specified in a network with missing dyads. Missingness will be overridden.")
        nw[as.matrix(is.na(nw),matrix.type="edgelist")] <- 0
      }else if(obs.constraints!=~.-observed){
        cat("Target statistics specified in a network with a nontrivial observation process. Observation process will be ignored.\n")
      }
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
