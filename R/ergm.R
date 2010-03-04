ergm <- function(formula, theta0="MPLE",
                 MPLEonly=FALSE, MLestimate=!MPLEonly, seed=NULL,
                 burnin=10000, MCMCsamplesize=10000, interval=100,
                 maxit=3,
                 constraints=~.,
                 meanstats=NULL,
                 dissolve=NULL, gamma=-4.59512, dissolve.order="FormAndDiss", # this line not in CRAN
                 control=control.ergm(),
                 verbose=FALSE, ...) {
  current.warn <- options()$warn
  options(warn=0)
  if(!is.null(seed))  set.seed(as.integer(seed))
  if (verbose) cat("Evaluating network in model\n")

  nw <- ergm.getnetwork(formula)
  if(!is.null(meanstats)){
   control$drop <- FALSE
   if(!(!is.null(control$SAN.burnin) && is.na(control$SAN.burnin))){
    netsumm<-summary(formula)
    if(length(netsumm)!=length(meanstats))
      stop("Incorrect length of the meanstats vector: should be ", length(netsumm), " but is ",length(meanstats),".")
    
    if(verbose) cat("Constructing an approximate response network.\n")
    ## If meanstats are given, overwrite the given network and formula
    ## with SAN-ed network and formula.
    nw<-san(formula, meanstats=meanstats,
            theta0=if(is.numeric(theta0)) theta0, 
            constraints=constraints,
            verbose=verbose,
            burnin=
            if(is.null(control$SAN.burnin)) burnin
            else control$SAN.burnin,
            interval=interval)
    formula<-ergm.update.formula(formula,nw~.)
    if (verbose) {
     cat("Original meanstats:\n")
     print(meanstats)
     cat("SAN meanstats - Original meanstats:\n")
     print(summary(formula, basis=nw)-meanstats)
    }
   }
  }
  if(control$nsubphases=="maxit") control$nsubphases<-maxit
  
  if (verbose) cat("Initializing model.\n")

  proposalclass <- if(is.null(dissolve)) "c" else "f"  # This line not in CRAN version
    
  if(control$drop){
   model.initial <- ergm.getmodel(formula, nw, drop=FALSE, initialfit=TRUE)
   model.initial.drop <- ergm.getmodel(formula, nw, drop=TRUE, initialfit=TRUE)
   namesmatch <- match(model.initial$coef.names, model.initial.drop$coef.names)
   droppedterms <- rep(FALSE, length=length(model.initial$etamap$offsettheta))
   droppedterms[is.na(namesmatch)] <- TRUE
   model.initial$etamap$offsettheta[is.na(namesmatch)] <- TRUE
  }else{
   model.initial <- ergm.getmodel(formula, nw, drop=control$drop, initialfit=TRUE)
   droppedterms <- rep(FALSE, length=length(model.initial$etamap$offsettheta))
  }
  if (verbose) cat("Initializing Metropolis-Hastings proposal.\n")
  MHproposal <- MHproposal(constraints, weights=control$prop.weights, control$prop.args, nw, model.initial,class=proposalclass)
  # Note:  MHproposal function in CRAN version does not use the "class" argument for now
  MHproposal.miss <- MHproposal("randomtoggleNonObserved", control$prop.args, nw, model.initial)

  conddeg <- switch(MHproposal$name=="CondDegree",control$drop,NULL)
  MCMCparams=c(control,
   list(samplesize=MCMCsamplesize, burnin=burnin, interval=interval,
        maxit=maxit,Clist.miss=NULL, mcmc.precision=control$mcmc.precision))


  if (verbose) cat("Fitting initial model.\n")
  theta0copy <- theta0
  initialfit <- ergm.initialfit(theta0=theta0copy, MLestimate=MLestimate, 
                                formula=formula, nw=nw, meanstats=meanstats,
				m=model.initial,
                                MPLEtype=control$MPLEtype, 
                                initial.loglik=control$initial.loglik,
                                conddeg=conddeg, MCMCparams=MCMCparams, MHproposal=MHproposal,
                                force.MPLE=(ergm.independencemodel(model.initial)
                                            && constraints==(~.)),
				verbose=verbose, 
                                compressflag = control$compress, 
                                maxNumDyadTypes=control$maxNumDyadTypes,
                                ...)
  MCMCflag <- ((MLestimate && (!ergm.independencemodel(model.initial)
                               || !is.null(meanstats)
                               || constraints!=(~.)))
                || control$force.mcmc || !is.null(dissolve))
  if (MCMCflag) {
    theta0 <- initialfit$coef
    names(theta0) <- model.initial$coef.names
    theta0[is.na(theta0)] <- 0
  } else { # Just return initial (non-MLE) fit and exit.
    initialfit$offset <- model.initial$etamap$offsettheta
    initialfit$drop <- droppedterms
    initialfit$network <- nw
    initialfit$newnetwork <- nw
    initialfit$formula <- formula
    initialfit$constraints <- constraints
    initialfit$prop.args <- control$prop.args
    initialfit$prop.weights <- control$prop.weights
    return(initialfit)
  } 
  if(control$drop){
   model <- ergm.getmodel(formula, nw, drop=FALSE, expanded=TRUE, silent=TRUE)
   # revise theta0 to reflect additional parameters
   theta0 <- ergm.revisetheta0(model, theta0)
   model.drop <- ergm.getmodel(formula, nw, drop=TRUE, expanded=TRUE, silent=TRUE)
#            silent="MPLE" %in% theta0copy)
   namesdrop <- model$coef.names[is.na(match(model$coef.names, model.drop$coef.names))]
   names(model$etamap$offsettheta) <- names(theta0)
   droppedterms <- rep(FALSE, length=length(model$etamap$offsettheta))
   droppedterms[is.na(namesmatch)] <- TRUE
   theta0[droppedterms] <- -Inf
   model$etamap$offsettheta[names(model$etamap$offsettheta) %in% namesdrop] <- TRUE
  }else{
   model <- ergm.getmodel(formula, nw, drop=control$drop, expanded=TRUE)
   # revise theta0 to reflect additional parameters
   theta0 <- ergm.revisetheta0(model, theta0)
   names(model$etamap$offsettheta) <- names(theta0)
   droppedterms <- rep(FALSE, length=length(model$etamap$offsettheta))
  }

  Clist <- ergm.Cprepare(nw, model)
  Clist.miss <- ergm.design(nw, model, verbose=verbose)
  Clist$obs <- summary(model$formula, drop=FALSE)
# Clist$obs <- summary(model$formula, drop=control$drop)
  Clist$meanstats <- Clist$obs
  if(!is.null(meanstats)){
   if (is.null(names(meanstats))){
    if(length(meanstats) == length(Clist$obs)){
     names(meanstats) <- names(Clist$obs)
     Clist$meanstats <- meanstats
    }else{
     namesmatch <- names(summary(model$formula, drop=FALSE))
     if(length(meanstats) == length(namesmatch)){
       namesmatch <- match(names(meanstats), namesmatch)
       Clist$meanstats <- meanstats[namesmatch]
     }
    }
   }else{
    namesmatch <- match(names(Clist$obs), names(meanstats))
    Clist$meanstats[!is.na(namesmatch)] <- meanstats[namesmatch[!is.na(namesmatch)]]
   }
  }

  MCMCparams=c(control,
   list(samplesize=MCMCsamplesize, burnin=burnin, interval=interval,
        maxit=maxit,Clist.miss=Clist.miss, mcmc.precision=control$mcmc.precision))

  if(!is.null(dissolve)){  # This section not in CRAN version.
    if (verbose) cat("Fitting Dynamic ERGM.\n")
    dissolve<-ergm.update.formula(dissolve,nw~.)
    model.dissolve <- ergm.getmodel(dissolve, nw, dissolve.order=dissolve.order)
    MHproposal.diss <- MHproposal(constraints, weights=control$prop.weights.diss, control$prop.args.diss, nw, model.dissolve,class="d")
    v <- switch(control$style.dyn,
                "SPSA" = ergm.SPSA.dyn(theta0, nw, model, model.dissolve,
                  Clist, gamma, 
                  MCMCparams=MCMCparams, MHproposal.form=MHproposal,
                  MHproposal.diss=MHproposal.diss,MT=FALSE,
                  verbose),
                "SPSA2" = ergm.SPSA.dyn(theta0, nw, model, model.dissolve,
                  Clist, gamma, 
                  MCMCparams=MCMCparams, MHproposal.form=MHproposal,
                  MHproposal.diss=MHproposal.diss,MT=TRUE,
                  verbose),
                "Robbins-Monro" = ergm.robmon.dyn(theta0, nw, model, model.dissolve,
                  Clist, gamma, 
                  MCMCparams=MCMCparams, MHproposal.form=MHproposal,
                  MHproposal.diss=MHproposal.diss,
                  verbose),
                ergm.mainfitloop.dyn(theta0, nw,
                                     model.form=model, model.diss=model.dissolve, Clist,
                                     gamma, initialfit,
                                     MCMCparams=MCMCparams, 
                                     MHproposal.form=MHproposal, MHproposal.diss=MHproposal.diss,
                                     verbose=verbose, 
                                     ...)
                )
  }else{
   if (verbose) cat("Fitting ERGM.\n")
   v <- switch(control$style,
    "Robbins-Monro" = ergm.robmon(theta0, nw, model, Clist, burnin, interval,
                      MHproposal(constraints,weights=control$prop.weights, control$prop.args, nw, model), verbose, control),
    "PILA" = ergm.PILA(theta0, nw, model, Clist,
      MHproposal(constraints,weights=control$prop.weights, control$prop.args, nw, model),
      MCMCparams, control, verbose),
    #  PILA stuff:  Not in CRAN version
    "Stochastic-Approximation" = ergm.stocapprox(theta0, nw, model, 
                                 Clist, 
                                 MCMCparams=MCMCparams, MHproposal=MHproposal,
                                 verbose),
                      ergm.mainfitloop(theta0, nw,
                          model, Clist, 
                          initialfit,
                          MCMCparams=MCMCparams, MHproposal=MHproposal,
                          MHproposal.miss=MHproposal.miss,
                          verbose=verbose, 
                          ...)
              )
  }
  if(!is.null(MCMCparams$check.degeneracy) && MCMCparams$check.degeneracy && (is.null(v$theta1$independent) || !all(v$theta1$independent))){
    if(verbose) {
      cat("Checking for degeneracy.\n")
    }
    degeneracy <- ergm.degeneracy(v, test.only=TRUE)
  } else {
    degeneracy <- list(degeneracy.value=NULL, degeneracy.type=NULL)
  }
  v$degeneracy.value <- degeneracy$degeneracy.value
  v$degeneracy.type <- degeneracy$degeneracy.type
  v$formula <- formula
  v$formula.diss <- dissolve  # Not in CRAN version
  v$constraints <- constraints
  v$prop.args <- control$prop.args
  v$prop.weights <- control$prop.weights
  v$prop.args.diss <- control$prop.args.diss # Not in CRAN version
  v$prop.weights.diss <- control$prop.weights.diss # Not in CRAN version

  v$offset <- model$etamap$offsettheta
  v$drop <- droppedterms
  v$etamap <- model$etamap
  if (!control$returnMCMCstats)
    v$sample <- NULL
  options(warn=current.warn)
  if (MCMCflag) {
    cat("\nThis model was fit using MCMC.  To examine model diagnostics", 
        "and check for degeneracy, use the mcmc.diagnostics() function.\n")
  }
  v
}
