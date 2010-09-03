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
  proposalclass <- "c"
  # Next for conditional MLE in dynamic model 
  if(is.character(MLestimate) && 
     (MLestimate=="formation") | (MLestimate=="dissolution")){
   lhs <- terms(formula)[[2]]
   if( is.call(lhs) && (lhs[[1]]=="|")){
   if(MLestimate=="formation"){
    proposalclass <- "fmle"
    y0 <- ergm.getnetwork(as.formula(paste("~",lhs[[3]])))
    y1 <- ergm.getnetwork(as.formula(paste("~",lhs[[2]])))
    nw = network.copy(y0)
    nwm = network.copy(y0)
    set.edge.attribute(nwm,attrname="na",value=TRUE)
    ydiscordantedges <- as.matrix(network(as.sociomatrix(y1)&!as.sociomatrix(y0),directed=is.directed(y0)),matrix.type="edgelist")
    add.edges(nw,ydiscordantedges[,1],ydiscordantedges[,2])
    add.edges(nwm,ydiscordantedges[,1],ydiscordantedges[,2])
#   delete.edges(y0,eid=1:length(y0$mel))
#   add.edges(y0,ydiscordantedges[,1],ydiscordantedges[,2])
   }
   if(MLestimate=="dissolution"){
    proposalclass <- "dmle"
    y0 <- ergm.getnetwork(as.formula(paste("~",lhs[[3]])))
    y1 <- ergm.getnetwork(as.formula(paste("~",lhs[[2]])))
#   nw = network.copy(y0)
#   delete.edges(nw,eid=1:length(nw$mel))
    nw = copy.null.network(y0)
    nwm <- as.sociomatrix(y0)
    y0nonedges <- as.matrix(network((nwm==0)&col(nwm)!=row(nwm),directed=is.directed(y0)),matrix.type="edgelist")
    nwm = copy.null.network(y0)
    add.edges(nwm,y0nonedges[,1],y0nonedges[,2])
    set.edge.attribute(nwm,attrname="na",value=TRUE)
    yminusedges <- as.matrix(network(as.sociomatrix(y0)&as.sociomatrix(y1),directed=is.directed(y0)),matrix.type="edgelist")
    add.edges(nw,yminusedges[,1],yminusedges[,2])
    add.edges(nwm,yminusedges[,1],yminusedges[,2])
#   ydiscordantedges <- as.matrix(network(as.sociomatrix(y0)&!as.sociomatrix(y1),directed=is.directed(y0)),matrix.type="edgelist")
#   delete.edges(y0,eid=1:length(y0$mel))
#   add.edges(y0,ydiscordantedges[,1],ydiscordantedges[,2])
   }
   formula.passed<-formula
   formula<-ergm.update.formula(formula,nw~.)
   MLestimate=!MPLEonly
  }}else{
   nwm <- nw
  }
  # End conditional MLE in dynamic model
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
                                formula=formula, nw=nwm, meanstats=meanstats,
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
  if((MHproposal$name!="FormationMLE")&(MHproposal$name!="DissolutionMLE")){
    Clist.miss <- ergm.design(nw, model, verbose=verbose)
  }else{
    Clist.miss <- ergm.Cprepare(y0, model)
  }
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
    "PILA" = ergm.PILA(theta0, nw, model, Clist,
      MHproposal(constraints,weights=control$prop.weights, control$prop.args, nw, model),
      MCMCparams, control, verbose),
    #  PILA stuff:  Not in CRAN version
    "Stochastic-Approximation" = ergm.stocapprox(theta0, nw, model, 
                                 Clist, 
                                 MCMCparams=MCMCparams, MHproposal=MHproposal,
                                 verbose),
    "Stepping" = ergm.stepping(theta0, nw, model, Clist, initialfit, 
				#nstats=nstats, 
				#approx=lognormapprox, filename.prefix=NULL, 
				#control=control.ergm(nsim1=100, nsim2=1000, gridsize=100),  # simulation parameters
				#plots=FALSE,  # currently useless, but plots can be reimplemented
				MCMCparams=MCMCparams, 
				MHproposal=MHproposal, MHproposal.miss=MHproposal.miss, 
				verbose=verbose,...),
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
  if(exists("formula.passed")){
    v$formula <- formula.passed
  }else{
    v$formula <- formula
  }
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
