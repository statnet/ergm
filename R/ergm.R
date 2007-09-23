ergm <- function(formula, theta0="MPLE", 
                 MPLEonly=FALSE, MLestimate=!MPLEonly, seed=NULL,
                 burnin=10000, MCMCsamplesize=10000, interval=100,
                 maxit=3,
                 constraint="none",
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
    
  model.initial <- ergm.getmodel(formula, nw, drop=control$drop, initialfit=TRUE)
  MHproposal <- getMHproposal(lookupMHproposal(proposalclass,constraint,control$prop.weights), control$prop.args, nw, model.initial)
  MHproposal.miss <- getMHproposal("randomtoggleNonObserved", control$prop.args, nw, model.initial)

  # MPLE & Meanstats -> need fake network
  if("MPLE" %in% theta0 && !is.null(meanstats)){
    nw.initial<-mk.pseudonet(meanstats, formula, nw, verbose=verbose)
  }
  else nw.initial<-nw
  
  BD <- ergm.boundDeg(control$boundDeg, nnodes=network.size(nw))
  Clist.initial <- ergm.Cprepare(nw.initial, model.initial)
  Clist.miss.initial <- ergm.design(nw.initial, model.initial, initialfit=TRUE,
                                verbose=verbose)
  Clist.initial$meanstats=meanstats
  theta0copy <- theta0
  initialfit <- ergm.initialfit(theta0copy, MLestimate, Clist.initial,
                                Clist.miss.initial, model.initial, verbose=verbose, ...)
  if (MLestimate && 
      (   !ergm.independencemodel(model.initial)
       || !is.null(meanstats))
       || control$force.mcmc
      ) {
    theta0 <- initialfit$coef
    names(theta0) <- model.initial$coef.names
    theta0[is.na(theta0)] <- 0
  } else { # Just return initial (non-MLE) fit and exit.
    initialfit$network <- nw
    initialfit$newnetwork <- nw
    initialfit$formula <- formula
    return(initialfit)
  } 
  model <- ergm.getmodel(formula, nw, drop=control$drop, expanded=TRUE)
  theta0 <- ergm.revisetheta0(model, theta0)
  # revise theta0 to reflect additional parameters

  Clist <- ergm.Cprepare(nw, model)
  Clist.miss <- ergm.design(nw, model, verbose=verbose)
  Clist$obs <- summary(model$formula, drop=control$drop)
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
  MCMCparams=c(control,list(samplesize=MCMCsamplesize, burnin=burnin, interval=interval,maxit=maxit,Clist.miss=Clist.miss))
  styles <- c("Newton-Raphson","Robbins-Monro","Stochastic-Approximation")
  control$style <- styles[pmatch(control$style,styles,nomatch=1)]
  if(!is.null(dissolve)){
    if (verbose) cat("Fitting Dynamic ERGM.\n")
    model.dissolve <- ergm.getmodel.dissolve(dissolve, nw, dissolve.order)
    MHproposal.diss <- getMHproposal(lookupMHproposal("d",constraint,control$weights.diss), control$prop.args.diss, nw, model.dissolve)
    v <- switch(control$style,
                "Robbins-Monro" = ergm.robmon.dyn(theta0, nw, model, model.dissolve,
                  Clist, BD, gamma, 
                  MCMCparams=MCMCparams, MHproposal.form=MHproposal,
                  MHproposal.diss=MHproposal.diss,
                  verbose),
                ergm.mainfitloop.dyn(theta0, nw,
                                     model.form=model, model.diss=model.dissolve, Clist,
                                     BD, gamma, initialfit,
                                     MCMCparams=MCMCparams, 
                                     MHproposal.form=MHproposal, MHproposal.diss=MHproposal.diss,
                                     verbose=verbose, 
                                     ...)
                )
  }else{
   if (verbose) cat("Fitting ERGM.\n")
   v <- switch(control$style,
    "Robbins-Monro" = ergm.robmon(theta0, nw, model, Clist, BD, burnin, interval,
                      getMHproposal(lookupMHproposal("c",constraint,control$prop.weights), control$prop.args, nw, model), verbose, control),
    "Stochastic-Approximation" = ergm.stocapprox(theta0, nw, model, 
                                 Clist, BD, 
                                 MCMCparams=MCMCparams, MHproposal=MHproposal,
                                 verbose),
                      ergm.mainfitloop(theta0, nw,
                          model, Clist, Clist.miss,
                          BD, initialfit,
                          MCMCparams=MCMCparams, MHproposal=MHproposal,
                          MHproposal.miss=MHproposal.miss,
                          verbose=verbose, 
                          ...)
              )
  }

  if(!is.null(v$mplefit)){
   if(v$loglikelihood>control$trustregion-0.0001){
    v$degeneracy <- control$trustregion
   }else{
     ## Workaround.
     #fff <- (-2*v$mplefit$glm$y+1)*model.matrix(v$mplefit$glm)
     #v$degeneracy.type <- apply(fff,1, ergm.degeneracy,theta0, model, v$sample)
     #v$degeneracy <- max(v$degeneracy.type)
     v$degeneracy<-NA

   }
  }else{
   if(MPLEonly){
    v$degeneracy.type <- abs(model.matrix(fit$glm) %*% fit$coef)
    v$degeneracy <- max(v$degeneracy.type)
   }else{
    v$degeneracy <- control$trustregion
   }
  }

  if (!control$returnMCMCstats)
    v$sample <- NULL
  options(warn=current.warn)
  v
}
