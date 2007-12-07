ergm <- ergm2 <- ergm.ihs <- function(formula, theta0="MPLE", 
                 MPLEonly=FALSE, MLestimate=!MPLEonly, seed=NULL,
                 burnin=10000, MCMCsamplesize=10000, interval=100,
                 maxit=3,
                 constraints=~.,
                 meanstats=NULL,
                 dissolve=NULL, gamma=-4.59512, dissolve.order="DissThenForm",
                 control=ergm.control(),
                 verbose=FALSE, ...) {
  current.warn <- options()$warn
  options(warn=0)
  if(!is.null(seed))  set.seed(as.integer(seed))
  if (verbose) cat("Evaluating network in model\n")

  nw <- ergm.getnetwork(formula)
  if(!is.null(meanstats)){ control$drop <- FALSE }
  if(control$nsubphases=="maxit") control$nsubphases<-maxit

  
  if (verbose) cat("Fitting initial model.\n")

  proposalclass <- if(is.null(dissolve)) "c" else "f"
    
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
  MHproposal <- MHproposal(constraints, weights=control$prop.weights, control$prop.args, nw, model.initial,class=proposalclass)
  MHproposal.miss <- MHproposal("randomtoggleNonObserved", control$prop.args, nw, model.initial)

  # MPLE & Meanstats -> need fake network
  if("MPLE" %in% theta0 && !is.null(meanstats)){
  # if IHS 
    nw.initial<-san(formula, meanstats=meanstats, verbose=verbose)
  #else
  # nw.initial<-mk.pseudonet(meanstats, formula, nw, verbose=verbose)
  # IHS end
  }
  else nw.initial<-nw
  
  Clist.initial <- ergm.Cprepare(nw.initial, model.initial)
  Clist.miss.initial <- ergm.design(nw.initial, model.initial, initialfit=TRUE,
                                verbose=verbose)
  Clist.initial$meanstats=meanstats
  theta0copy <- theta0
  initialfit <- ergm.initialfit(theta0copy, MLestimate, Clist.initial,
                                Clist.miss.initial, model.initial,
                                MPLEtype=control$MPLEtype, 
                                initial.loglik=control$initial.loglik,
                                verbose=verbose, ...)
  if (MLestimate && 
      (   !ergm.independencemodel(model.initial)
       || !is.null(meanstats)
       || constraints!=(~.))
      || control$force.mcmc
      ) {
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

  if (verbose) cat("ergm.mainfitloop\n")
  MCMCparams=c(control,
   list(samplesize=MCMCsamplesize, burnin=burnin, interval=interval,maxit=maxit,Clist.miss=Clist.miss))
  if(!is.null(dissolve)){
    if (verbose) cat("Fitting Dynamic ERGM.\n")
    model.dissolve <- ergm.getmodel.dissolve(dissolve, nw, dissolve.order)
    MHproposal.diss <- MHproposal(constraints, weights=control$prop.weights.diss, control$prop.args.diss, nw, model.dissolve,class="d")
    v <- switch(control$style,
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
    "Stochastic-Approximation" = ergm.stocapprox(theta0, nw, model, 
                                 Clist, 
                                 MCMCparams=MCMCparams, MHproposal=MHproposal,
                                 verbose),
                      ergm.mainfitloop(theta0, nw,
                          model, Clist, Clist.miss,
                          initialfit,
                          MCMCparams=MCMCparams, MHproposal=MHproposal,
                          MHproposal.miss=MHproposal.miss,
                          verbose=verbose, 
                          ...)
              )
  }

  v$formula <- formula
  v$formula.diss <- dissolve
  v$constraints <- constraints
  v$prop.args <- control$prop.args
  v$prop.weights <- control$prop.weights
  v$prop.args.diss <- control$prop.args.diss
  v$prop.weights.diss <- control$prop.weights.diss

  v$offset <- model$etamap$offsettheta
  v$drop <- droppedterms
  v$etamap <- model$etamap
  if (!control$returnMCMCstats)
    v$sample <- NULL
  options(warn=current.warn)
  v
}
