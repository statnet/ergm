#  File R/ergm.model.utils.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################

# A helper function to check for extreme statistics and inform the user.
ergm.checkextreme.model <- function(model, nw, init, response, target.stats, drop, silent=FALSE){
  eta0 <- ergm.eta(init, model$etamap)
  
  # Check if any terms are at their extremes.
  obs.stats.eta <- if(!is.null(target.stats)){
    if(is.null(names(target.stats))) names(target.stats) <- names(ergm.getglobalstats(nw, model, response=response))
    target.stats
  }else ergm.getglobalstats(nw, model, response=response)

  extremeval.eta <- +(model$maxval==obs.stats.eta)-(model$minval==obs.stats.eta)
  names.eta<-names(obs.stats.eta)      
  
  # Drop only works for canonical terms, so extremeval.theta has 0s
  # (no action) for all elements of theta that go into a curved term,
  # and only the elements of extremeval.eta corresponding to
  # canonical terms get copied into it.
  extremeval.theta <- rep(0, length(init))
  extremeval.theta[model$etamap$canonical!=0 & !model$etamap$offsettheta]<-extremeval.eta[model$etamap$canonical[!model$etamap$offsettheta]]
  names.theta <- rep(NA, length(length(init)))
  names.theta[model$etamap$canonical!=0]<-names.eta[model$etamap$canonical]
  
  # The offsettheta is there to prevent dropping of offset terms.
  low.drop.theta <- extremeval.theta<0 & !model$etamap$offsettheta
  high.drop.theta <- extremeval.theta>0 & !model$etamap$offsettheta
  
  # drop only affects the action taken.
  if(drop){

    if(!silent){
      # Inform the user what's getting dropped.
      if(any(low.drop.theta)) message(paste("Observed statistic(s)", paste.and(names.theta[low.drop.theta]), "are at their smallest attainable values. Their coefficients will be fixed at -Inf.", sep=" "))
      if(any(high.drop.theta)) message(paste("Observed statistic(s)", paste.and(names.theta[high.drop.theta]), "are at their greatest attainable values. Their coefficients will be fixed at +Inf.", sep=" "))
    }
    
    # If the user specified a non-fixed element of init, and that element is getting dropped, warn the user.
    if(any(is.finite(init[low.drop.theta|high.drop.theta]))) warning("Overriding user-specified initial init coefficient", paste.and(names.theta[is.na(init[low.drop.theta|high.drop.theta])]), ". To preserve, enclose in an offset() function.", sep="")

    init[low.drop.theta|high.drop.theta] <- extremeval.theta[low.drop.theta|high.drop.theta]*Inf
    model$etamap$offsettheta[low.drop.theta|high.drop.theta] <- TRUE
    model$etamap$offsetmap[model$etamap$canonical[low.drop.theta|high.drop.theta & (model$etamap$canonical>0)]] <- TRUE
  }else{
    if(!silent){
      # If no drop, warn the user anyway.
      if(any(low.drop.theta)) warning(paste("Observed statistic(s)", paste.and(names.theta[low.drop.theta]), "are at their smallest attainable values and drop=FALSE. The MLE is poorly defined.", sep=" "))
      if(any(high.drop.theta)) warning(paste("Observed statistic(s)", paste.and(names.theta[high.drop.theta]), "are at their greatest attainable values and drop=FALSE. The MLE is poorly defined.", sep=" "))
    }
  }
  
  list(model=model, init=init, extremeval.theta=extremeval.theta)
}

## Check for conflicts between model terms and constraints.

ergm.checkconstraints.model <- function(model, MHproposal, init, silent=FALSE){
  # Get the list of all the constraints that the proposal imposes on the sample space.
  constraints.old<-names(MHproposal$arguments$constraints)
  repeat{
    constraints <- unique(sort(c(constraints.old, unlist(ergm.ConstraintImplications()[constraints.old]))))
    if(all(constraints %in% constraints.old)) break
    else constraints.old <- constraints
  }

  coef.counts <- coef.sublength.model(model)
  conflict.coefs <- c()
  
  for(i in seq_along(model$terms))
    conflict.coefs <- c(conflict.coefs, rep(any(model$terms[[i]]$conflicts.constraints %in% constraints), # Is there a conflict?
                                            coef.counts[i])) # How many coefficients are affected?

  conflict.coefs[model$offsettheta] <- FALSE # No conflict if it's already an offset.
  
  if(any(conflict.coefs)){
    if(!silent){
      warning(paste("The specified model's sample space constraint holds statistic(s)", paste.and(model$coef.names[conflict.coefs]), " constant. They will be ignored.", sep=" "))
    }
    init[conflict.coefs] <- 0
    model$etamap$offsettheta[conflict.coefs] <- TRUE
    model$etamap$offsetmap[model$etamap$canonical[conflict.coefs & (model$etamap$canonical>0)]] <- TRUE    
  }
  
  list(model=model, init=init, estimable=!conflict.coefs)
}

coef.sublength.model<-function(object, ...){
  sapply(object$terms, function(term){
    ## curved term
    if(!is.null(term$params)) length(term$params)
    ## linear term
    else length(term$coef.names)
  })
}

coef.length.model <- function(object, ...){
  sum(coef.sublength.model(object))
}

.coef.names.model <- function(object, canonical){
    if(canonical) object$coef.names
    else unlist(lapply(object$terms, function(term) NVL(names(term$params),term$coef.names)))
}
