ergm <- function(formula, theta0="MPLE", 
                 MPLEonly=FALSE, MLestimate=!MPLEonly, seed=NULL,
                 burnin=10000, MCMCsamplesize=10000, interval=100,
                 maxit=3,
                 constraints=~.,
                 control=control.ergm(),      
                 verbose=FALSE, ...) {
  current.warn <- options()$warn
  options(warn=0)
  if(!is.null(seed))  set.seed(as.integer(seed))
  if (verbose) cat("Evaluating network in model\n")
  nw <- ergm.getnetwork(formula)
  if(control$nsubphases=="maxit") control$nsubphases<-maxit  
                              
  
  if (verbose) cat("Initializing model.\n")    
   model.initial <- ergm.getmodel(formula, nw, drop=FALSE, initialfit=TRUE)
   droppedterms <- rep(FALSE, length=length(model.initial$etamap$offsettheta))
  if(control$drop){
   model.initial.drop <- ergm.getmodel(formula, nw, drop=TRUE, initialfit=TRUE)
   namesmatch <- match(model.initial$coef.names, model.initial.drop$coef.names)
   droppedterms[is.na(namesmatch)] <- TRUE
   model.initial$etamap$offsettheta[is.na(namesmatch)] <- TRUE
  }

  if (verbose) cat("Initializing Metropolis-Hastings proposal.\n")
  MHproposal <- MHproposal(constraints, weights=control$prop.weights, control$prop.args, nw, model.initial)

  if(!is.null(control$initial.network)){
    nw.initial<-control$initial.network
  }else{
    nw.initial<-nw
  }
  Clist.initial <- ergm.Cprepare(nw.initial, model.initial)
  Clist2.initial <- list(heads=0, tails=0, nedges=0, dir=is.directed(nw)) #unused for now

  if (verbose) cat("Fitting initial model.\n")
  Clist.initial$meanstats=NULL
  theta0copy <- theta0
  initialfit <- ergm.initialfit(theta0copy, MLestimate, Clist.initial,
                                Clist2.initial, model.initial,
                                MPLEtype=control$MPLEtype, verbose=verbose, ...)
  if (MLestimate && (!ergm.independencemodel(model.initial)
                     || constraints!=(~.))
     || control$force.mcmc) {
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
  Clist2 <- list(heads=0, tails=0, nedges=0, dir=is.directed(nw)) #unused for now
  Clist$obs <- summary(model$formula, drop=FALSE)
# Clist$obs <- summary(model$formula, drop=control$drop)
  Clist$meanstats <- Clist$obs

  MCMCparams=c(control,
   list(samplesize=MCMCsamplesize, burnin=burnin, interval=interval,maxit=maxit,Clist.miss=Clist2))
  if (verbose) cat("Fitting ERGM.\n")
    v <- switch(control$style,
                "Robbins-Monro" = ergm.robmon(theta0, nw, model, Clist, burnin, interval,
                                              MHproposal(constraints,
                                                            weights=control$prop.weights, 
                                                            control$prop.args, nw, model), 
                                              verbose, control),
                "Stochastic-Approximation" = ergm.stocapprox(theta0, nw, model, 
                                                             Clist, 
                                                             MCMCparams=MCMCparams, MHproposal=MHproposal,
                                                             verbose),
                ergm.mainfitloop(theta0, nw,
                                 model, Clist,
                                 initialfit,
                                 MCMCparams=MCMCparams, MHproposal=MHproposal,
                                 verbose=verbose, 
                                 ...)
                )
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
  v
}

