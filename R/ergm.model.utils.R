# A helper function to reorder vector v (if named) into order specified by names.
vector.namesmatch<-function(v,names,errname=NULL){
  if(is.null(errname)) errname <- deparse(substitute(v))

  if (is.null(names(v))){
    if(length(v) == length(names)){
      names(v) <- names
    }else stop("Length of ``", errname, "'' is ", length(v), " should be ", names,".")
  }else{
    if(length(v) == length(names)
       && length(unique(names(v)))==length(v)
       && length(unique(names))==length(names)
       && all(sort(names(v)) == sort(names))){
      namesmatch <- match(names(v), names)
      v <- v[namesmatch]
    }else stop("Name missmatch in ``", errname,"''. Specify by position.")
  }
  v
}

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
  extremeval.theta[model$etamap$canonical!=0]<-extremeval.eta[model$etamap$canonical]
  names.theta <- rep(NA, length(length(init)))
  names.theta[model$etamap$canonical!=0]<-names.eta[model$etamap$canonical]
  
  # The offsettheta is there to prevent dropping of offset terms.
  low.drop.theta <- extremeval.theta<0 & !model$etamap$offsettheta
  high.drop.theta <- extremeval.theta>0 & !model$etamap$offsettheta
  
  # drop only affects the action taken.
  if(drop){

    if(!silent){
      # Inform the user what's getting dropped.
      if(any(low.drop.theta)) cat(paste("Observed statistic(s)", paste(names.theta[low.drop.theta],collapse=", "), "are at their smallest attainable values. Their coefficients will be fixed at -Inf.\n", sep=" "))
      if(any(high.drop.theta)) cat(paste("Observed statistic(s)", paste(names.theta[high.drop.theta],collapse=", "), "are at their greatest attainable values. Their coefficients will be fixed at +Inf.\n", sep=" "))
    }
    
    # If the user specified a non-fixed element of init, and that element is getting dropped, warn the user.
    if(any(is.finite(init[low.drop.theta|high.drop.theta]))) warning("Overriding user-specified initial init coefficient", paste(names.theta[is.na(init[low.drop.theta|high.drop.theta])], collapse=", "), ". To preserve, enclose in an offset() function.", sep="")

    init[low.drop.theta|high.drop.theta] <- extremeval.theta[low.drop.theta|high.drop.theta]*Inf
    model$etamap$offsettheta[low.drop.theta|high.drop.theta] <- TRUE
    model$etamap$offsetmap[model$etamap$canonical[low.drop.theta|high.drop.theta & (model$etamap$canonical>0)]] <- TRUE
  }else{
    if(!silent){
      # If no drop, warn the user anyway.
      if(any(low.drop.theta)) warning(paste("Observed statistic(s)", paste(names.theta[low.drop.theta],collapse=", "), "are at their smallest attainable values and drop=FALSE. The MLE is poorly defined.", sep=" "))
      if(any(high.drop.theta)) warning(paste("Observed statistic(s)", paste(names.theta[high.drop.theta],collapse=", "), "are at their greatest attainable values and drop=FALSE. The MLE is poorly defined.", sep=" "))
    }
  }

  list(model=model, init=init, extremeval.theta=extremeval.theta)
}

## Check for conflicts between model terms and constraints.

ergm.checkconstraints.model <- function(model, MHproposal){
  # Get the list of all the constraints that the proposal imposes on the sample space.
  constraints.old<-names(MHproposal$arguments$constraints)
  repeat{
    constraints <- unique(sort(c(constraints.old, unlist(ConstraintImplications[constraints.old]))))
    if(all(constraints==constraints.old)) break
    else constraints.old <- constraints
  }

  conflict.coefs <- c()
  
  for(term in model$terms){
    if(any(term$conflicts.constraints %in% constraints))
      conflict.coefs <- c(conflict.coefs, term$coef.names)
  }

  if(length(conflict.coefs))
    warning(paste("The specified model's sample space constraint holds statistic(s)", paste(conflict.coefs, collapse=", "), " constant. They will be ignored.", sep=" "))
}

coef.sublength.model<-function(m){
  sapply(m$terms, function(term){
    ## curved term
    if(!is.null(term$params)) length(term$params)
    ## linear term
    else length(term$coef.names)
  })
}

coef.length.model<-function(m){
  sum(coef.sublength.model(m))
}
