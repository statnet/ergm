#  File R/ergm_model.utils.R in package ergm, part of the Statnet suite
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

#' @rdname coef.length.model
#'
#' @description \code{coef.sublength.model} returns a vector
#'   containing the number of model parameters corresponding to each
#'   model term.
#' @export coef.sublength.model
coef.sublength.model<-function(object, offset=NA, ...){
  terms <-
    if(is.na(offset)) object$terms
    else if(offset) object$terms[object$etamap$offset]
    else if(!offset) object$terms[!object$etamap$offset]
                  
  sapply(terms, function(term){
    ## curved term
    if(!is.null(term$params)) length(term$params)
    ## linear term
    else length(term$coef.names)
  })
}

#' Extract Number of parameters in ergm Model
#' 
#' \code{coef.length.model} returns the total number of parameters of
#' an [`ergm_model`].
#' 
#' @param object an [`ergm_model`] object
#' @param offset If `NA` (the default), all model terms are counted;
#'   if \code{TRUE}, only offset terms are counted; and if
#'   \code{FALSE}, offset terms are skipped.
#' @param \dots other arguments.
#' @return \code{coef.length.model} returns the length of the
#'   parameter vector of the model.
#' @note These are *not* methods at this time. This may change in the
#'   future.
#' @keywords models
#' @export coef.length.model
coef.length.model <- function(object, offset=NA, ...){
  sum(coef.sublength.model(object, offset=offset, ...))
}

#' @rdname coef.length.model
#' @description \code{eta.sublength.model} returns a vector containing
#'   the number of canonical parameters (also the number of sufficient
#'   statistics) corresponding to each model term.
#' @export eta.sublength.model
eta.sublength.model<-function(object, offset=NA, ...){
  terms <-
    if(is.na(offset)) object$terms
    else if(offset) object$terms[object$etamap$offset]
    else if(!offset) object$terms[!object$etamap$offset]
                  
  sapply(terms, function(term){
    length(term$coef.names)
  })
}

#' @rdname coef.length.model
#' @description \code{eta.length.model} returns the total number of
#'   canonical parameters.
#' @export eta.sublength.model
eta.length.model <- function(object, offset=NA, ...){
  sum(eta.sublength.model(object, offset=offset, ...))
}

#' Parameters names of an initialized ergm model
#'
#' @param object an `ergm_model`.
#' @param canonical whether the canonical (natural) parameters or the
#'   model (curved) parameters are wanted.
#'
#' @return a character vector of parameter names.
coef.names.model <- function(object, canonical){
  if(canonical) object$coef.names
  else unlist(lapply(object$terms, function(term) NVL(names(term$params),term$coef.names)))
}

#' `ergm_model`'s `etamap` with all offset terms removed and remapped
#'
#' @param etamap the `etamap` element of the `ergm_model`.
#'
#' @return a copy of the input as if all offset terms were dropped
#'   from the model
deoffset.etamap <- function(etamap){
  if(!any(etamap$offsetmap)) return(etamap) # Nothing to do.

  # Construct "mappings" from model with offset to model without.
  remap.theta <- cumsum(!etamap$offsettheta)
  remap.theta[etamap$offsettheta] <- 0
  remap.eta <- cumsum(!etamap$offsetmap)
  remap.eta[etamap$offsetmap] <- 0

  out <- list()  
  # Canonical map: for those non-offset elements of $canonical that
  # are not 0, copy them over after remapping them to their new
  # positions.
  out$canonical <- rep(0, sum(!etamap$offsettheta))
  out$canonical[etamap$canonical[!etamap$offsettheta]>0] <- remap.eta[etamap$canonical[!etamap$offsettheta]]

  # These are just the non-FALSE elements.
  out$offsetmap <- etamap$offsetmap[!etamap$offsetmap]
out$offsettheta <- etamap$offsettheta[!etamap$offsettheta]
  out$offset <- etamap$offset[!etamap$offset]

  # For each curved term,
  out$curved <- lapply(etamap$curved, function(term){
    term$from <- remap.theta[term$from] # remap from
    term$to <- remap.eta[term$to] # remap to
    term
  })
  # Filter out those terms for which the "to" map is not all nonzeros.
  out$curved <- out$curved[!sapply(lapply(lapply(out$curved, "[[", "to"),"==",0), any)]
  out$etalength <- length(out$offsetmap)

  out
}