#  File R/obs.constraints.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################

.handle.obs.constraints <- function(nw,
                                    constraints=~.,
                                    obs.constraints=~observed,
                                    target.stats=NULL) {
  
  # Observation process handling only needs to happen if the
  # sufficient statistics are not specified. If the sufficient
  # statistics are specified, the nw's dyad states are irrelevant.
  if(!is.null(target.stats)){
    if(network.naedgecount(nw)){
      warning("Target statistics specified in a network with missing dyads. Missingness will be overridden.")
      nw[as.matrix(is.na(nw),matrix.type="edgelist")] <- 0
    }else if(obs.constraints!=~observed){
      cat("Target statistics specified in a network with a nontrivial observation process. Observation process will be ignored.\n")
    }
    obs.constraints <- ~.
  }
  # Get list of observation process constraints, dropping the . in ~. if needed.
  obs.constraints <- term.list.formula(obs.constraints[[length(obs.constraints)]])
  obs.constraints <- obs.constraints[obs.constraints!="."]

  # If no missing edges, remove the "observed" constraint.
  if(network.naedgecount(nw)==0){
    obs.con.names <- sapply(obs.constraints, function(x) as.character(if(is.call(x)) x[[1]] else x))
    obs.constraints[obs.con.names=="observed"] <- NULL
  }
  
  constraints.obs<-append.rhs.formula(constraints, obs.constraints, TRUE)
  if(constraints==constraints.obs) constraints.obs<-NULL

  list(nw = nw, constraints.obs = constraints.obs)
}

.align.target.stats.offset <- function(model, target.stats){
  om <- model$etamap$offsetmap
  cno <- coef.names.mode(model, canonical=TRUE)
  cn <- coef.names.mode(model, canonical=TRUE)[!om]
  target.stats <- na.omit(vector.namesmatch(target.stats, cn))
  tmp <- rep(NA, length(om))
  tmp[!om] <- target.stats
  names(tmp) <- cno
  tmp
}
