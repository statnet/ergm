#  File ergm/R/ergm.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of Washington
#                David R. Hunter, Penn State University
#                Carter T. Butts, University of California - Irvine
#                Steven M. Goodreau, University of Washington
#                Martina Morris, University of Washington
# Copyright 2007 The statnet Development Team
######################################################################
ergm <- function(formula, theta0="MPLE",
                 MPLEonly=FALSE, MLestimate=!MPLEonly, seed=NULL,
                 burnin=10000, MCMCsamplesize=10000, interval=100,
                 maxit=3,
                 constraints=~.,
                 meanstats=NULL,
                 control=control.ergm(),
                 verbose=FALSE, ...) {
  current.warn <- options()$warn
  options(warn=0)
  if(!is.null(seed))  set.seed(as.integer(seed))
  if (verbose) cat("Evaluating network in model\n")

  nw <- ergm.getnetwork(formula)
  if(!is.null(meanstats)){
    netsumm<-summary(formula)
    if(length(netsumm)!=length(meanstats))
      stop("Incorrect length of the meanstats vector: should be ", length(netsumm), " but is ",length(meanstats),".")
    
    control$drop <- FALSE

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
    formula<-safeupdate.formula(formula,nw~.)
    if (verbose) {
     cat("Original meanstats:\n")
     print(meanstats)
     cat("Original meanstats - SAN meanstats:\n")
     print(summary(formula, basis=nw)-meanstats)
    }
  }
  if(control$nsubphases=="maxit") control$nsubphases<-maxit
  
  if (verbose) cat("Initializing model.\n")

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
  MHproposal <- MHproposal(constraints, weights=control$prop.weights, control$prop.args, nw, model.initial)
  MHproposal.miss <- MHproposal("randomtoggleNonObserved", control$prop.args, nw, model.initial)

  Clist.initial <- ergm.Cprepare(nw, model.initial)
  Clist.miss.initial <- ergm.design(nw, model.initial, verbose=verbose)
  if (verbose) cat("Fitting initial model.\n")
  Clist.initial$meanstats=meanstats
  theta0copy <- theta0
  initialfit <- ergm.initialfit(theta0copy, MLestimate, Clist.initial,
                                Clist.miss.initial, model.initial,
                                MPLEtype=control$MPLEtype, 
                                initial.loglik=control$initial.loglik,
                                verbose=verbose, compressflag = control$compress, 
                                maxNumDyadTypes=control$maxNumDyadTypes,
                                force.MPLE=(ergm.independencemodel(model.initial)
                                            && constraints==(~.)),
                                ...)
  MCMCflag <- ((MLestimate && (!ergm.independencemodel(model.initial)
                               || !is.null(meanstats)
                               || constraints!=(~.)))
                || control$force.mcmc)
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

   if (verbose) cat("Fitting ERGM.\n")
   v <- switch(control$style,
    "Robbins-Monro" = ergm.robmon(theta0, nw, model, Clist, burnin, interval,
                      MHproposal(constraints,weights=control$prop.weights, control$prop.args, nw, model), verbose, control),
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
  v$constraints <- constraints
  v$prop.args <- control$prop.args
  v$prop.weights <- control$prop.weights

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
