stergm <- function(formation, dissolution, theta.form0="MPLE", theta.diss=NULL,
                   seed=NULL,
                   MH.burnin=1000,
                   constraints=~., # May not work.
                   meanstats=NULL,
                   stergm.order="FormAndDiss",
                 control=control.stergm(),
                 verbose=FALSE, ...) {
  current.warn <- options()$warn
  options(warn=0)
  if(!is.null(seed))  set.seed(as.integer(seed))
  if (verbose) cat("Evaluating network in model.form\n")

  nw <- ergm.getnetwork(formation)
  if(!is.null(meanstats)){
   control$drop <- FALSE
   netsumm<-summary(formation)
   if(length(netsumm)!=length(meanstats))
     stop("Incorrect length of the meanstats vector: should be ", length(netsumm), " but is ",length(meanstats),".")
   
   if(verbose) cat("Constructing an approximate response network.\n")
   ## If meanstats are given, overwrite the given network and formation
   ## with SAN-ed network and formation.
   nw<-san(formation, meanstats=meanstats,
           theta0=if(is.numeric(theta.form0)) theta.form0, 
           constraints=constraints,
           verbose=verbose,
           burnin=control$SAN.burnin,
           interval=control$SAN.interval)
   formation<-ergm.update.formula(formation,nw~.)
   if (verbose) {
     cat("Original meanstats:\n")
     print(meanstats)
     cat("SAN meanstats - Original meanstats:\n")
     print(summary(formation, basis=nw)-meanstats)
   }
  }
  
  if (verbose) cat("Initializing model.\n")

  # Note: I am not sure whether this works:
  if(control$drop){
   model.initial <- ergm.getmodel(formation, nw, drop=FALSE, initialfit=TRUE)
   model.initial.drop <- ergm.getmodel(formation, nw, drop=TRUE, initialfit=TRUE)
   namesmatch <- match(model.initial$coef.names, model.initial.drop$coef.names)
   droppedterms <- rep(FALSE, length=length(model.initial$etamap$offsettheta))
   droppedterms[is.na(namesmatch)] <- TRUE
   model.initial$etamap$offsettheta[is.na(namesmatch)] <- TRUE
  }else{
   model.initial <- ergm.getmodel(formation, nw, drop=control$drop, initialfit=TRUE)
   droppedterms <- rep(FALSE, length=length(model.initial$etamap$offsettheta))
  }
  
  if (verbose) cat("Initializing Metropolis-Hastings proposal.\n")
  MHproposal.form <- MHproposal(constraints, weights=control$prop.weights.form, control$prop.args.form, nw, model.initial,class="f")

  conddeg <- switch(MHproposal.form$name=="CondDegree",control$drop,NULL)
  MCMCparams=c(control,
   list(MH.burnin=MH.burnin))


  if (verbose) cat("Fitting initial model.\n")
  theta0copy <- theta.form0
  initialfit <- ergm.initialfit(theta0=theta0copy, MLestimate=TRUE, 
                                formula=formation, nw=nw, meanstats=meanstats,
                                m=model.initial,
                                MPLEtype=control$MPLEtype, 
                                initial.loglik=control$initial.loglik,
                                conddeg=conddeg, MCMCparams=MCMCparams, MHproposal=MHproposal.form,
                                force.MPLE=FALSE,
                                verbose=verbose, 
                                compressflag = control$compress, 
                                maxNumDyadTypes=control$maxNumDyadTypes,
                                ...)

  theta.form0 <- initialfit$coef
  names(theta.form0) <- model.initial$coef.names
  theta.form0[is.na(theta.form0)] <- 0

  
  if(control$drop){
   model.form <- ergm.getmodel(formation, nw, drop=FALSE, expanded=TRUE, silent=TRUE)
   # revise theta.form0 to reflect additional parameters
   theta.form0 <- ergm.revisetheta0(model.form, theta.form0)
   model.drop <- ergm.getmodel(formation, nw, drop=TRUE, expanded=TRUE, silent=TRUE)
#            silent="MPLE" %in% theta.form0copy)
   namesdrop <- model.form$coef.names[is.na(match(model.form$coef.names, model.drop$coef.names))]
   names(model.form$etamap$offsettheta) <- names(theta.form0)
   droppedterms <- rep(FALSE, length=length(model.form$etamap$offsettheta))
   droppedterms[is.na(namesmatch)] <- TRUE
   theta.form0[droppedterms] <- -Inf
   model.form$etamap$offsettheta[names(model.form$etamap$offsettheta) %in% namesdrop] <- TRUE
  }else{
   model.form <- ergm.getmodel(formation, nw, drop=control$drop, expanded=TRUE)
   # revise theta.form0 to reflect additional parameters
   theta.form0 <- ergm.revisetheta0(model.form, theta.form0)
   names(model.form$etamap$offsettheta) <- names(theta.form0)
   droppedterms <- rep(FALSE, length=length(model.form$etamap$offsettheta))
  }

  Clist <- ergm.Cprepare(nw, model.form)
  Clist$obs <- summary(model.form$formula, drop=FALSE)
  Clist$meanstats <- Clist$obs
  if(!is.null(meanstats)){
   if (is.null(names(meanstats))){
    if(length(meanstats) == length(Clist$obs)){
     names(meanstats) <- names(Clist$obs)
     Clist$meanstats <- meanstats
    }else{
     namesmatch <- names(summary(model.form$formula, drop=FALSE))
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
   list(MH.burnin=MH.burnin))

if (verbose) cat("Fitting Dynamic ERGM.\n")
    dissolution<-ergm.update.formula(dissolution,nw~.)
    model.diss <- ergm.getmodel(dissolution, nw, stergm.order=stergm.order)
    MHproposal.diss <- MHproposal(constraints, weights=control$prop.weights.diss, control$prop.args.diss, nw, model.diss,class="d")
    v <- switch(control$style,
                "SPSA" = stergm.SPSA(theta.form0, nw, model.form, model.diss,
                  Clist, theta.diss, 
                  MCMCparams=MCMCparams, MHproposal.form=MHproposal.form,
                  MHproposal.diss=MHproposal.diss,MT=FALSE,
                  verbose),
                "SPSA2" = stergm.SPSA(theta.form0, nw, model.form, model.diss,
                  Clist, theta.diss, 
                  MCMCparams=MCMCparams, MHproposal.form=MHproposal.form,
                  MHproposal.diss=MHproposal.diss,MT=TRUE,
                  verbose),
                "Nelder-Mead" = stergm.NM(theta.form0, nw, model.form, model.diss,
                  Clist, theta.diss, 
                  MCMCparams=MCMCparams, MHproposal.form=MHproposal.form,
                  MHproposal.diss=MHproposal.diss,
                  verbose),
                "Robbins-Monro" = stergm.RM(theta.form0, nw, model.form, model.diss,
                  Clist, theta.diss, 
                  MCMCparams=MCMCparams, MHproposal.form=MHproposal.form,
                  MHproposal.diss=MHproposal.diss,
                  verbose),
                stergm.mainfitloop(theta.form0, nw,
                                     model.form=model.form, model.diss=model.diss, Clist,
                                     theta.diss, initialfit,
                                     MCMCparams=MCMCparams, 
                                     MHproposal.form=MHproposal.form, MHproposal.diss=MHproposal.diss,
                                     verbose=verbose, 
                                     ...)
                )

  v$MH.burnin <- MH.burnin
  v$stergm.order <- stergm.order
  v$formation <- formation
  v$dissolution <- dissolution
  v$constraints <- constraints
  v$prop.args.form <- control$prop.args.form
  v$prop.weights.form <- control$prop.weights.form
  v$prop.args.diss <- control$prop.args.diss
  v$prop.weights.diss <- control$prop.weights.diss

  v$offset <- model.form$etamap$offsettheta
  v$drop <- droppedterms
  v$etamap <- model.form$etamap
  options(warn=current.warn)
  v
}
