#  File R/ergm_model.utils.R in package ergm, part of the Statnet suite of
#  packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################

# A helper function to check for extreme statistics and inform the user.
ergm.checkextreme.model <- function(model, nw, init, response, target.stats, drop, silent=FALSE){
  eta0 <- ergm.eta(init, model$etamap)
  
  # Check if any terms are at their extremes.
  obs.stats.eta <- if(!is.null(target.stats)){
    if(is.null(names(target.stats))) names(target.stats) <- param_names(model, canonical=TRUE)
    target.stats
  }else summary(model, nw)

  extremeval.eta <- +(model$maxval==obs.stats.eta)-(model$minval==obs.stats.eta)
  names.eta<-names(obs.stats.eta)      
  
  # Drop only works for canonical terms, so extremeval.theta has 0s
  # (no action) for all elements of theta that go into a curved term,
  # and only the elements of extremeval.eta corresponding to
  # canonical terms get copied into it.
  extremeval.theta <- rep_along(init, 0)
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

ergm.checkconstraints.model <- function(model, proposal, init, silent=FALSE){
  # Get the list of all the constraints that the proposal imposes on the sample space.
  # FIXME: This doesn't cover "second-order" constraints. Is it worth
  # it to bring back the implications API?
  constraints <- unique(sort(as.character(unlist(lapply(proposal$arguments$constraints, `[[`, "implies")))))
  
  coef.counts <- nparam(model, byterm=TRUE)
  conflict.coefs <- c()
  
  for(i in seq_along(model$terms))
    conflict.coefs <- c(conflict.coefs, rep(any(model$terms[[i]]$conflicts.constraints %in% constraints), # Is there a conflict?
                                            coef.counts[i])) # How many coefficients are affected?

  conflict.coefs[model$offsettheta] <- FALSE # No conflict if it's already an offset.
  
  if(any(conflict.coefs)){
    if(!silent){
      warning(paste("The specified model's sample space constraint holds statistic(s)", paste.and(param_names(model,canonical=TRUE)[conflict.coefs]), " constant. They will be ignored.", sep=" "))
    }
    init[conflict.coefs] <- 0
    model$etamap$offsettheta[conflict.coefs] <- TRUE
    model$etamap$offsetmap[model$etamap$canonical[conflict.coefs & (model$etamap$canonical>0)]] <- TRUE    
  }
  
  list(model=model, init=init, estimable=!conflict.coefs)
}

#' `ergm_model`'s `etamap` with all offset terms removed and remapped
#'
#' @param etamap the `etamap` element of the `ergm_model`.
#' @param coef the vector of model coefficients specifying the offsets, if any
#'
#' @return a copy of the input as if all offset terms were dropped
#'   from the model; if a term is a partial offset, the term is kept
#'   but the map and gradient functions are rewritten to substitute
#'   the elements of `coef`. It should be safe to call with
#'   `theta=coef[!etamap$offsettheta]`.
#' @noRd
deoffset.etamap <- function(etamap, coef){
  if(!any(etamap$offsettheta)) return(etamap) # Nothing to do.

  # Construct "mappings" from model with offset to model without.
  remap.theta <- as.integer(cumsum(!etamap$offsettheta))
  remap.theta[etamap$offsettheta] <- 0L
  remap.eta <- as.integer(cumsum(!etamap$offsetmap))
  remap.eta[etamap$offsetmap] <- 0L

  out <- list()  
  # Canonical map: for those non-offset elements of $canonical that
  # are not 0, copy them over after remapping them to their new
  # positions.
  out$canonical <- integer(sum(!etamap$offsettheta))
  out$canonical[etamap$canonical[!etamap$offsettheta]>0L] <- remap.eta[etamap$canonical[!etamap$offsettheta]]

  # For each curved term,
  out$curved <- lapply(etamap$curved, function(term){
    offtheta <- etamap$offsettheta[term$from]
    if(all(offtheta)) return(NULL) # Drop the whole term

    if(any(offtheta)){ # Only some are offsets...
      mymap <- term$map
      mygrad <- term$gradient
      myx <- coef[term$from]

      term$map <- function(x, ...){
        myx[!offtheta] <- x
        mymap(myx, ...)
      }
      term$gradient <- function(x, ...){
        myx[!offtheta] <- x
        mygrad(myx, ...)[!offtheta,,drop=FALSE]
      }

      term$from <- term$from[!offtheta]
    }

    term$from <- remap.theta[term$from] # remap from
    term$to <- remap.eta[term$to] # remap to
    term
  }) %>% compact()

  # These are just the non-FALSE elements.
  out$offsetmap <- etamap$offsetmap[!etamap$offsetmap]
  out$offsettheta <- etamap$offsettheta[!etamap$offsettheta]

  out$etalength <- length(out$offsetmap)

  out
}
