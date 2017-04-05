#  File ergm/R/stergm.EGMME.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
################################################################################
# The <stergm> function fits stergms from a specified formation and dissolution
# formula returning approximate MLE's based on MCMC estimation.
################################################################################

stergm.EGMME <- function(nw, formation, dissolution,  offset.coef.form, offset.coef.diss,
                   targets, target.stats, estimate,
                 control,
                 verbose) {

  if(!is.network(nw)) stop("Argument nw must be a network.")
   
  # Allow the user to specify targets as copied from formation or dissolution formula.
  if(is.character(targets)){
    targets <- switch(targets,
                      formation = formation,
                      dissolution = dissolution)
  }
  
  if(length(targets)==3){
    warning("Targets formula has an LHS, which will be ignored in favor of nw.")
    targets <- targets[c(1,3)] # in a formula f<-y~x, f[1]=~, f[2]=y, and f[3]=x
  }

  targets <- ergm.update.formula(targets,nw~.)
  formation <- ergm.update.formula(formation,nw~.)
  dissolution <- ergm.update.formula(dissolution,nw~.)
  
  if(!is.null(target.stats)){
    nw.stats<-summary(targets)
    if(length(nw.stats)!=length(target.stats))
      stop("Incorrect length of the target.stats vector: should be ", length(nw.stats), " but is ",length(target.stats),".")
    
    if(verbose) cat("Constructing an approximate response network.\n")
    ## If target.stats are given, overwrite the given network and targets
    ## with SAN-ed network and targets.
    newnw <- try({
      for(srun in seq_len(control$SAN.maxit)){
        nw<-san(targets, target.stats=target.stats,
                control=control$SAN.control,
                verbose=verbose)
        targets<-ergm.update.formula(targets,nw~.)
        nw.stats <- summary(targets)
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
      nw
    }, silent=TRUE)
    if(inherits(newnw,"try-error")){
      cat("SAN failed or is not applicable. Increase burn-in if there are problems.\n")
    }else nw <- newnw
    formation <- ergm.update.formula(formation,nw~.)
    dissolution <- ergm.update.formula(dissolution,nw~.)
  }

  if (verbose) cat("Initializing Metropolis-Hastings proposals.\n")
  MHproposal.form <- MHproposal(~., weights=control$MCMC.prop.weights.form, control$MCMC.prop.args.form, nw, class="f")
  MHproposal.diss <- MHproposal(~., weights=control$MCMC.prop.weights.diss, control$MCMC.prop.args.diss, nw, class="d")
  
  model.form <- ergm.getmodel(formation, nw, expanded=TRUE, role="formation")
  model.diss <- ergm.getmodel(dissolution, nw, expanded=TRUE, role="dissolution")
  model.mon <- ergm.getmodel(targets, nw, expanded=TRUE, role="target")

  if(any(model.form$etamap$canonical==0) || any(model.diss$etamap$canonical==0) || any(model.mon$etamap$canonical==0)) stop("Equilibrium GMME for models based on curved ERGMs is not supported at this time.")

  p.free<-sum(!model.form$etamap$offsettheta)+sum(!model.diss$etamap$offsettheta)
  if(p.free==0) stop("Model specification has no free parameters (all are offsets).")
  q<-length(model.mon$etamap$offsettheta)
  if(p.free>q) stop("Fitting ",p.free," free parameters on ",q," target statistics. The specification is underidentified.")

  model.mon$nw.stats <- summary(model.mon$formula)
  model.mon$target.stats <- if(!is.null(target.stats)) vector.namesmatch(target.stats, names(model.mon$nw.stats)) else model.mon$nw.stats

  # If some control$init is specified...
  
  if(!is.null(control$init.form)){
    # Check length of control$init.form.
    if(length(control$init.form)!=length(model.form$etamap$offsettheta)) {
      if(verbose) cat("control$init.form is", control$init.form, "\n", "number of statistics is",length(model.form$coef.names), "\n")
      stop(paste("Invalid starting formation parameter vector control$init.form:",
                 "wrong number of parameters."))
    }
  }else control$init.form <- rep(NA, length(model.form$etamap$offsettheta)) # Set the default value of control$init.form.
  if(!is.null(offset.coef.form)) control$init.form[model.form$etamap$offsettheta]<-offset.coef.form
  names(control$init.form) <- model.form$coef.names

  if(!is.null(control$init.diss)){
    # Check length of control$init.diss.
    if(length(control$init.diss)!=length(model.diss$etamap$offsettheta)) {
      if(verbose) cat("control$init.diss is", control$init.diss, "\n", "number of statistics is",length(model.diss$coef.names), "\n")
      stop(paste("Invalid starting dissolution parameter vector control$init.diss:",
                 "wrong number of parameters."))
    }
  }else control$init.diss <- rep(NA, length(model.diss$etamap$offsettheta)) # Set the default value of control$init.diss.  
  if(!is.null(offset.coef.diss)) control$init.diss[model.diss$etamap$offsettheta]<-offset.coef.diss
  names(control$init.diss) <- model.diss$coef.names

  initialfit <- stergm.EGMME.initialfit(control$init.form, control$init.diss, nw, model.form, model.diss, model.mon, control, verbose)
  
  if(verbose) cat("Fitting STERGM Equilibrium GMME.\n")

  Cout <- switch(control$EGMME.main.method,
                 "Stochastic-Approximation" = stergm.EGMME.SA(initialfit$formation.fit$coef,
                   initialfit$dissolution.fit$coef, nw, model.form, model.diss, model.mon,
                   control=control, MHproposal.form=MHproposal.form,
                  MHproposal.diss=MHproposal.diss,
                  verbose),
                 stop("Method ", control$EGMME.main.method, " is not implemented.")
                )

  out <- list(network = nw, formation = formation, dissolution = dissolution, targets = targets, target.stats=model.mon$target.stats, estimate=estimate, covar = Cout$covar, opt.history=Cout$opt.history, sample=Cout$sample, sample.obs=NULL, control=control,
              formation.fit = with(Cout, list(network=nw, formula=formation, coef = eta.form, covar=covar.form, etamap = model.form$etamap, offset = model.form$etamap$offsettheta, constraints=~., estimate=estimate, control=control)),
              dissolution.fit = with(Cout, list(network=nw, formula=dissolution, coef = eta.diss, covar=covar.diss, etamap = model.diss$etamap, offset = model.diss$etamap$offsettheta, constraints=~., estimate=estimate, control=control))
              )
  class(out$formation.fit)<-class(out$dissolution.fit)<-"ergm"
  
  out
}
