################################################################################
# The <stergm> function fits stergms from a specified formation and dissolution
# formula returning approximate MLE's based on MCMC estimation.
#
# --PARAMETERS--
#   formation   : the formation formula, as 'nw ~ term(s)'
#   dissolution : the dissolution formula, as 'nw ~ term(s)'
#   theta.form0 : the intial theta formation coefficients, or optionally if
#                 these are to be estimates, the string "MPLE";
#                 default="MPLE"
#   theta.diss  : the initial theta dissolution coefficients
#   seed        : an integer starting value for the random number generator;
#                 default=NULL
#   MH.burnin   : the number of proposals used in each MCMC step; this is ignored
#                 unless 'control$main.method'="Robbins-Monro"; any other style or
#                 the default style will not recognize this parameter;
#                 default=1000
#   constraints : a one-sided formula of the constraint terms; options are
#                      bd        degrees        nodegrees
#                      edges     degreedist     indegreedist
#                      observed  outdegreedist
#                default="~ ."; these may not work currently.
#   target.stats   :  a vector of the mean value parameters;
#                  default=the observed statistic from the 'nw' in formula
#   control     :  a list of control parameters returned from <control.stergm>;
#                  default=<control.stergm>()
#   verbose     :  whether ergm should be verbose (T or F); default=FALSE
#
# --RETURNED--
#   because a stergm object is the return type of several functions, and
#   because this is a rather lengthy list, and because the returned items
#   of this function borrow from the other stergm.* functions, this list 
#   provides the returned items for all funtions returning a stergm.
#   The symbol preceding each component indicates which function returns it,
#   but remember that, <stergm> will additionally return the items from
#   one of the other stergm functions as well:
#       <stergm>             = $
#       <stergm.RM>          = @
#       <stergm.SPSA>        = &
#       <stergm.NM>          = !
#
#   the components include:
#
#     @   coef.form   : the estimated formation model coefficients
#     @   coef.diss   : the estimated dissolution model coefficients
#      &  eta         : the estimated formation ?? coefficients
#   $     offset      : a logical vector whose ith entry tells whether the
#                        ith curved theta coeffient was offset/fixed
#   $     etamap      :  the list constituting the theta->eta mapping for the
#                        formation model; for details of its components,
#                        see <ergm.etamap>
#   $     MH.burnin   :  the number of proposals made in each MCMC step
#   $     formation   : the formation formula, as 'nw ~ term(s)'
#   $     dissolution : the dissolution formula, as 'nw ~ term(s)'
#   $     constraints : the constraints formula
#     @&  newnetwork  :  the final network sampled
#     @&  network    :  the 'nw' inputted to <ergm> via the 'formula'
#   $     prop.args.form     :  the MHP formation arguments passed to the
#                               InitMHP rountines
#   $     prop.args.diss     :  the MHP dissolution arguments passed to the
#                               InitMHP rountines
#   $     prop.weights.form  :  the method used to allocate probabilities of
#                               being proposed to dyads in the formation stage,
#                               as "TNT", "random", "nonobserved", or "default"
#      &  theta.original     :  the theta values at the start of the MCMC 
#                               sampling
#     @   theta.form.original:  the formation theta values at the start of the
#                               MCMC sampling
#   $     prop.weights.diss  :  as 'prop.weights.form', but for the dissolution
#                               model
#     &   objective.history  :  the number of SPSA iteration used
#
################################################################################

stergm.EGMoME <- function(nw, formation, dissolution,  offset.coef.form, offset.coef.diss,
                   targets, target.stats, estimate,
                 control,
                 verbose) {

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

  # FIXME: Offsets in formation are not handled.
  if(formation!=targets) stop("Fitting equillibrium GMoME with targets different from formation terms is not supported at this time.")

  targets <- ergm.update.formula(targets,nw~.)
  formation <- ergm.update.formula(formation,nw~.)
  dissolution <- ergm.update.formula(dissolution,nw~.)
  
  if(!is.null(target.stats)){
    netsumm<-summary(targets)
    if(length(netsumm)!=length(target.stats))
      stop("Incorrect length of the target.stats vector: should be ", length(netsumm), " but is ",length(target.stats),".")
    
    if(verbose) cat("Constructing an approximate response network.\n")
    ## If target.stats are given, overwrite the given network and targets
    ## with SAN-ed network and targets.
    srun <- 0
    obs <- target.stats-target.stats
    while((srun < control$SAN.maxit) & (sum((obs-target.stats)^2) > 5)){
      nw<-san(targets, target.stats=target.stats,
              constraints=~.,
              control=control$SAN.control,
              verbose=verbose)
      targets<-ergm.update.formula(targets,nw~.)
      obs <- summary(targets, basis=nw)
      srun <- srun + 1
      if(verbose){
        cat(paste("Finished SAN run",srun,"\n"))
      }
      if(verbose){
        cat("SAN summary statistics:\n")
        print(obs)
        cat("Meanstats Goal:\n")
        print(target.stats)
        cat("Difference: SAN target.stats - Goal target.stats =\n")
        print(round(obs-target.stats,0))
      }
    }
  }

  if (verbose) cat("Initializing Metropolis-Hastings proposals.\n")
  MHproposal.form <- MHproposal(~., weights=control$MCMC.prop.weights.form, control$MCMC.prop.args.form, nw, class="f")
  MHproposal.diss <- MHproposal(~., weights=control$MCMC.prop.weights.diss, control$MCMC.prop.args.diss, nw, class="d")
  
  model.form <- ergm.getmodel(formation, nw, expanded=TRUE, MHp=MHproposal.form)
  model.diss <- ergm.getmodel(dissolution, nw, expanded=TRUE, MHp=MHproposal.diss)

  if(any(model.form$etamap$canonical==0) && any(model.form$etamap$canonical==0)) stop("Equilibrium GMoME for curved ERGMs is not supported at this time.")
  if(!all(model.diss$offset)) stop("Equilibrium GMoME to estimate dissolution parameters is not supported at this time. All dissolution parameters must be offset().")

  model.form$obs <- summary(model.form$formula)
  model.form$target.stats <- if(!is.null(target.stats)) vector.namesmatch(target.stats, names(model.form$obs)) else model.form$obs

  # If some control$init is specified...
  
  if(!is.null(control$init.form)){
    # Check length of control$init.form.
    if (length(control$init.form)!=length(model.form$etamap$offsettheta)) {
      if(verbose) cat("control$init.form is", control$init, "\n", "number of statistics is",length(model.initial$coef.names), "\n")
      stop(paste("Invalid starting formation parameter vector control$init.form:",
                 "wrong number of parameters."))
    }
  }else control$init.form <- rep(NA, length(model.form$etamap$offsettheta)) # Set the default value of control$init.form.
  if(!is.null(offset.coef.form)) control$init.form[model.form$etamap$offsettheta]<-offset.coef.form
  names(control$init.form) <- model.form$coef.names

  if(!is.null(control$init.diss)){
    # Check length of control$init.diss.
    if (length(control$init.diss)!=length(model.diss$etamap$offsettheta)) {
      if(verbose) cat("control$init.diss is", control$init, "\n", "number of statistics is",length(model.initial$coef.names), "\n")
      stop(paste("Invalid starting dissolution parameter vector control$init.diss:",
                 "wrong number of parameters."))
    }
  }else control$init.diss <- rep(NA, length(model.diss$etamap$offsettheta)) # Set the default value of control$init.diss.  
  if(!is.null(offset.coef.diss)) control$init.diss[model.diss$etamap$offsettheta]<-offset.coef.diss
  names(control$init.diss) <- model.diss$coef.names

  initialfit <- stergm.EGMoME.initialfit(control$init.form, control$init.diss, nw, formation, dissolution, targets, target.stats, model.form, model.diss, control, verbose)
  
  if(verbose) cat("Fitting Dynamic ERGM.\n")

  Cout <- switch(control$EGMoME.main.method,
                "SPSA" = stergm.SPSA(initialfit$formation.fit$coef, nw, model.form, model.diss,
                  initialfit$dissolution.fit$coef, 
                  control=control, MHproposal.form=MHproposal.form,
                  MHproposal.diss=MHproposal.diss,MT=FALSE,
                  verbose),
                "SPSA2" = stergm.SPSA(initialfit$formation.fit$coef, nw, model.form, model.diss,
                  initialfit$dissolution.fit$coef, 
                  control=control, MHproposal.form=MHproposal.form,
                  MHproposal.diss=MHproposal.diss,MT=TRUE,
                  verbose),
                "Robbins-Monro" = stergm.RM(initialfit$formation.fit$coef, nw, model.form, model.diss,
                  initialfit$dissolution.fit$coef, 
                  control=control, MHproposal.form=MHproposal.form,
                  MHproposal.diss=MHproposal.diss,
                  verbose)
                )

  out <- list(network = nw, formation = formation, dissolution = dissolution, targets = targets, target.stats=target.stats, estimate=estimate,
              formation.fit = with(Cout, list(coef=coef.form, sample=NULL, sample.obs=NULL, objective.history=Cout$objective.history, etamap = model.form$etamap, offset = model.form$etamap$offsettheta)),
              dissolution.fit = list(coef = initialfit$dissolution.fit$coef, sample=NULL, sample.obs=NULL, etamap = model.diss$etamap, offset = model.diss$etamap$offsettheta))
  class(out$formation.fit)<-class(out$dissolution.fit)<-"ergm"
  
  out
}
