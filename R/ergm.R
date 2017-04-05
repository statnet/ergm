#  File ergm/R/ergm.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
###############################################################################
# The <ergm> function fits ergms from a specified formula returning either
# MPLEs or approximate MLE's based on MCMC estimation.
#####################################################################################    

ergm <- function(formula,
                 constraints=~.,
                 offset.coef=NULL,
                 target.stats=NULL,
                 eval.loglik=TRUE,
                 estimate=c("MLE", "MPLE"),
                 control=control.ergm(),
                 verbose=FALSE,...) {
  current.warn <- options()$warn
  on.exit(options(warn=current.warn), add=TRUE)
  options(warn=0)

  estimate <- match.arg(estimate)
  # Backwards-compatibility:
  control<-control.ergm.toplevel(control,...)
  if(!is.null(list(...)$MPLEonly) && list(...)$MPLEonly){
    warning("Argument MPLEonly is deprecated. Use ``estimate=\"MPLE\"'' instead." )
    estimate <- "MPLE"
  }

  
  if(!is.null(control$seed))  set.seed(as.integer(control$seed))
  if (verbose) cat("Evaluating network in model\n")
  
  nw <- ergm.getnetwork(formula)
  proposalclass <- "c"

  # Construct the constraint for the observation process.
  # There may be a better way to specify this in the future.
  
  MHproposal.obs<-constraints
  
  # Missing data handling only needs to happen if the sufficient
  # statistics are not specified. If the sufficient statistics are
  # specified, the nw's dyad states are irrelevant.
  if(network.naedgecount(nw)){
    if(!is.null(target.stats)){
      warning("Target statistics specified in a network with missing dyads. Missingness will be overridden.")
      nw[as.matrix(is.na(nw),matrix.type="edgelist")] <- 0
    }else MHproposal.obs<-ergm.update.formula(MHproposal.obs,~.+observed)
  }
  
  if(constraints==MHproposal.obs) MHproposal.obs<-NULL

  ## Construct approximate response network if target.stats are given.
  
  if(!is.null(target.stats)){
    nw.stats<-summary(formula)
    if(length(nw.stats)!=length(target.stats))
      stop("Incorrect length of the target.stats vector: should be ", length(nw.stats), " but is ",length(target.stats),".")
    
    if(verbose) cat("Constructing an approximate response network.\n")
    ## If target.stats are given, overwrite the given network and formula
    ## with SAN-ed network and formula.
    for(srun in 1:control$SAN.maxit){
      nw<-san(formula, target.stats=target.stats,
              constraints=constraints,
              control=control$SAN.control,
              verbose=verbose)
      formula<-ergm.update.formula(formula,nw~.)
      nw.stats <- summary(formula, basis=nw)
      srun <- srun + 1
      if(verbose){
        cat(paste("Finished SAN run",srun,"\n"))
      }
      if(verbose){
        cat("SAN summary statistics:\n")
        print(nw.stats)
        cat("Meanstats Goal:\n")
        print(target.stats)
        cat("Difference: SAN target.stats - Goal target.stats =\n")
        print(round(nw.stats-target.stats,0))
      }
      if(sum((nw.stats-target.stats)^2) <= 5) break
    }
  }
  
  if (verbose) { cat("Initializing Metropolis-Hastings proposal.\n") }
  
  MHproposal <- MHproposal(constraints, weights=control$MCMC.prop.weights, control$MCMC.prop.args, nw, class=proposalclass)
  # Note:  MHproposal function in CRAN version does not use the "class" argument for now
  if(!is.null(MHproposal.obs)) MHproposal.obs <- MHproposal(MHproposal.obs, weights=control$MCMC.prop.weights, control$MCMC.prop.args, nw, class=proposalclass)
  
  conddeg <- switch(MHproposal$name %in% c("CondDegree","CondDegreeSimpleTetrad","BipartiteCondDegHexadToggles","BipartiteCondDegTetradToggles"),control$drop,NULL)
  
  if (verbose) cat("Initializing model.\n")
  
  # Construct the initial model.
  model.initial <- ergm.getmodel(formula, nw, initialfit=TRUE)
   
  # If some control$init is specified...
  if(!is.null(control$init)){
    # Check length of control$init.
    if (length(control$init)!=length(model.initial$etamap$offsettheta)) {
      if(verbose) cat("control$init is", control$init, "\n", "number of statistics is",length(model.initial$coef.names), "\n")
      stop(paste("Invalid starting parameter vector control$init:",
                 "wrong number of parameters.",
                 "If you are passing output from another ergm run as control$init,",
                 "in a model with curved terms, see help(enformulate.curved)."))
    }
  }else control$init <- rep(NA, length(model.initial$etamap$offsettheta)) # Set the default value of control$init.

  if(!is.null(offset.coef)) control$init[model.initial$etamap$offsettheta]<-offset.coef
  
  # Make sure any offset elements are given in control$init.
  if(any(is.na(control$init) & model.initial$etamap$offsettheta)) stop("The model contains offset terms whose parameter values have not been specified:", paste.and(model.initial$coef.names[is.na(control$init)|model.initial$offsettheta]), ".", sep="")

  # Check if any terms are constrained to a constant and issue a warning.
  constrcheck <- ergm.checkconstraints.model(model.initial, MHproposal, control$init)
  model.initial <- constrcheck$model; control$init <- constrcheck$init

  # Check if any terms are at their extremes and handle them depending on control$drop.
  extremecheck <- ergm.checkextreme.model(model=model.initial, nw=nw, init=control$init, target.stats=target.stats, drop=control$drop)
  model.initial <- extremecheck$model; control$init <- extremecheck$init

  if (verbose) { cat("Fitting initial model.\n") }

  MPLE.is.MLE <- (ergm.independencemodel(model.initial)
                  && !control$force.main
                  && constraints==(~.))

  # If all other criteria for MPLE=MLE are met, _and_ SAN network matches target.stats directly, we can get away with MPLE.
  MCMCflag <- (estimate=="MLE" && (!MPLE.is.MLE
                               || (!is.null(target.stats) && !isTRUE(all.equal(target.stats,nw.stats)))
                              )
               || control$force.main)

  # Short-circuit the optimization if all terms are either offsets or dropped.
  if(MCMCflag && all(model.initial$etamap$offsettheta)){
    # Note that this cannot be overridden with control$force.main.
    MCMCflag <- FALSE
    warning("All terms are either offsets or extreme values. Skipping MCMC.")
  }
  
  initialfit <- ergm.initialfit(init=control$init, initial.is.final=!MCMCflag,
                                formula=formula, nw=nw, target.stats=target.stats,
                                m=model.initial, method=control$init.method,
                                MPLEtype=control$MPLE.type, 
                                conddeg=conddeg, control=control, MHproposal=MHproposal,
                                verbose=verbose, 
                                compressflag = control$MCMC.compress, 
                                maxNumDyadTypes=control$MPLE.max.dyad.types,
                                ...)
  
  


  if (MCMCflag) {
    init <- initialfit$coef
    names(init) <- model.initial$coef.names
    init[is.na(init)] <- 0
  } else { # Just return initial (non-MLE) fit and exit.
    initialfit$offset <- model.initial$etamap$offsettheta
    initialfit$drop <- if(control$drop) extremecheck$extremeval.theta
    initialfit$estimable <- constrcheck$estimable
    initialfit$network <- nw
    initialfit$newnetwork <- nw
    initialfit$formula <- formula
    initialfit$constrained <- MHproposal$arguments$constraints
    initialfit$constraints <- constraints
    initialfit$target.stats <- model.initial$target.stats
    initialfit$estimate <- estimate

    initialfit$control<-control
    
    if(any(!model.initial$etamap$offsettheta))
      initialfit<-logLik.ergm(initialfit, add=TRUE, control=control$loglik.control)
    return(initialfit)
  }
  
  # Construct the curved model
  model <- ergm.getmodel(formula, nw, expanded=TRUE, silent=TRUE)
  # revise init to reflect additional parameters
  init <- ergm.reviseinit(model, init)

  # Check if any terms are constrained to a constant and issue a warning.
  constrcheck <- ergm.checkconstraints.model(model, MHproposal, init=init, silent=TRUE)
  model <- constrcheck$model; control$init <- constrcheck$init
  
  # Check if any terms are at their extremes and handle them depending on control$drop.
  extremecheck <- ergm.checkextreme.model(model=model, nw=nw, init=init, target.stats=target.stats, drop=control$drop, silent=TRUE)
  model <- extremecheck$model; init <- extremecheck$init

  model$nw.stats <- summary(model$formula)
  model$target.stats <- if(!is.null(target.stats)) vector.namesmatch(target.stats, names(model$nw.stats)) else model$nw.stats

  if (verbose) cat("Fitting ERGM.\n")
  mainfit <- switch(control$main.method,
    "Robbins-Monro" = ergm.robmon(init, nw, model, 
                      MHproposal(constraints,weights=control$MCMC.prop.weights, control$MCMC.prop.args, nw), verbose, control),
    "Stochastic-Approximation" = ergm.stocapprox(init, nw, model, 
                                 control=control, MHproposal=MHproposal,
                                 verbose),
    "Stepping" = ergm.stepping(init, nw, model, initialfit, constraints,
				#nstats=nstats, 
				#approx=lognormapprox, filename.prefix=NULL, 
				#control=control.ergm(nsim1=100, nsim2=1000, gridsize=100),  # simulation parameters
				#plots=FALSE,  # currently useless, but plots can be reimplemented
				control=control, 
				MHproposal=MHproposal, MHproposal.obs=MHproposal.obs, 
				verbose=verbose,...),
    "MCMLE" = ergm.MCMLE(init, nw,
                          model, 
                          initialfit,
                          control=control, MHproposal=MHproposal,
                          MHproposal.obs=MHproposal.obs,
                          verbose=verbose,
                          ...),
              stop("Method ", control$main.method, " is not implemented.")
              )

  if(!is.null(control$MCMLE.check.degeneracy) && control$MCMLE.check.degeneracy && (is.null(mainfit$theta1$independent) || !all(mainfit$theta1$independent))){
    if(verbose) {
      cat("Checking for degeneracy.\n")
    }
    degeneracy <- ergm.degeneracy(mainfit, test.only=TRUE)
  } else {
    degeneracy <- list(degeneracy.value=NULL, degeneracy.type=NULL)
  }
  mainfit$degeneracy.value <- degeneracy$degeneracy.value
  mainfit$degeneracy.type <- degeneracy$degeneracy.type

  mainfit$formula <- formula
  mainfit$target.stats <- model$target.stats

  mainfit$constrained <- MHproposal$arguments$constraints
  mainfit$constraints <- constraints

  mainfit$control<-control

  mainfit$estimate <- estimate

  mainfit$offset <- model$etamap$offsettheta
  mainfit$drop <- if(control$drop) extremecheck$extremeval.theta
  mainfit$estimable <- constrcheck$estimable
  mainfit$etamap <- model$etamap
  if (!control$MCMC.return.stats)
    mainfit$sample <- NULL

  if (MCMCflag) {
    cat("\nThis model was fit using MCMC.  To examine model diagnostics", 
        "and check for degeneracy, use the mcmc.diagnostics() function.\n")
  }
  if(eval.loglik)
    mainfit<-logLik.ergm(mainfit, add=TRUE, control=control$loglik.control)
  mainfit
}
