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
#                 unless 'control$style'="Robbins-Monro"; any other style or
#                 the default style will not recognize this parameter;
#                 default=1000
#   constraints : a one-sided formula of the constraint terms; options are
#                      bd        degrees        nodegrees
#                      edges     degreedist     indegreedist
#                      observed  outdegreedist
#                default="~ ."; these may not work currently.
#   target.stats   :  a vector of the mean value parameters;
#                  default=the observed statistic from the 'nw' in formula
#   stergm.order:  the order in which the formation and dissolution processes
#                  should occur, as one of
#                      "DissThenForm"      "DissAndForm"
#                      "FormAndDiss"       "FormThenDiss"
#                      "FormOnly"           "DissOnly"
#                  default="FormAndDiss"
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
#   $     stergm.order:  the order in which the formation and dissolution processes
#                        were used, as one of
#                           "DissThenForm"      "DissAndForm"
#                           "FormAndDiss"       "FormThenDiss"
#                           "FormOnly"           "DissOnly"
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

stergm <- function(formation, dissolution, theta.form0=NULL, theta.diss=NULL,
                   seed=NULL,
                   MH.burnin=1000,
                   constraints=~., # May not work.
                   target.stats=NULL,
                   stergm.order="FormAndDiss",
                 control=control.stergm(),
                 verbose=FALSE, ...) {
  current.warn <- options()$warn
  options(warn=0)
  if(!is.null(seed))  set.seed(as.integer(seed))
  if (verbose) cat("Evaluating network in model.form\n")

  nw <- ergm.getnetwork(formation)
  if(!is.null(target.stats)){
   netsumm<-summary(formation)
   if(length(netsumm)!=length(target.stats))
     stop("Incorrect length of the target.stats vector: should be ", length(netsumm), " but is ",length(target.stats),".")
   
   if(verbose) cat("Constructing an approximate response network.\n")
   ## If target.stats are given, overwrite the given network and formation
   ## with SAN-ed network and formation.
   nw<-san(formation, target.stats=target.stats,
           theta0=if(is.numeric(theta.form0)) theta.form0, 
           constraints=constraints,
           verbose=verbose,
           burnin=control$SAN.burnin,
           interval=control$SAN.interval)
   formation<-ergm.update.formula(formation,nw~.)
   if (verbose) {
     cat("Original target.stats:\n")
     print(target.stats)
     cat("SAN target.stats - Original target.stats:\n")
     print(summary(formation, basis=nw)-target.stats)
   }
  }

  if (verbose) cat("Initializing Metropolis-Hastings proposal.\n")
  MHproposal.form <- MHproposal(constraints, weights=control$prop.weights.form, control$prop.args.form, nw, class="f")
  
  if (verbose) cat("Initializing model.\n")

  model.initial <- ergm.getmodel(formation, nw, initialfit=TRUE, MHp=MHproposal.form)
  if(is.null(theta.form0)) theta.form0 <- rep(NA, length(model.initial$etamap$offsettheta))
  

  MCMCparams=c(control,
   list(MH.burnin=MH.burnin))


  if (verbose) cat("Fitting initial model.\n")
  initialfit <- ergm.initialfit(theta0=theta.form0, initial.is.final=FALSE, 
                                formula=formation, nw=nw, target.stats=target.stats,
                                m=model.initial, method=control$initialfit,
                                MPLEtype=control$MPLEtype, 
                                MCMCparams=MCMCparams, MHproposal=MHproposal.form,
                                verbose=verbose, 
                                compressflag = control$compress, 
                                maxNumDyadTypes=control$maxNumDyadTypes,
                                ...)

  theta.form0 <- initialfit$coef
  names(theta.form0) <- model.initial$coef.names
  theta.form0[is.na(theta.form0)] <- 0

  
   model.form <- ergm.getmodel(formation, nw, expanded=TRUE, MHp=MHproposal.form)
   # revise theta.form0 to reflect additional parameters
   theta.form0 <- ergm.revisetheta0(model.form, theta.form0)
   names(model.form$etamap$offsettheta) <- names(theta.form0)

  Clist <- ergm.Cprepare(nw, model.form)
  Clist$obs <- summary(model.form$formula)
  Clist$target.stats <- Clist$obs
  if(!is.null(target.stats)){
   if (is.null(names(target.stats))){
    if(length(target.stats) == length(Clist$obs)){
     names(target.stats) <- names(Clist$obs)
     Clist$target.stats <- target.stats
    }else{
     namesmatch <- names(summary(model.form$formula))
     if(length(target.stats) == length(namesmatch)){
       namesmatch <- match(names(target.stats), namesmatch)
       Clist$target.stats <- target.stats[namesmatch]
     }
    }
   }else{
    namesmatch <- match(names(Clist$obs), names(target.stats))
    Clist$target.stats[!is.na(namesmatch)] <- target.stats[namesmatch[!is.na(namesmatch)]]
   }
  }

  MCMCparams=c(control,
   list(MH.burnin=MH.burnin))

if (verbose) cat("Fitting Dynamic ERGM.\n")
  dissolution<-ergm.update.formula(dissolution,nw~.)
  MHproposal.diss <- MHproposal(constraints, weights=control$prop.weights.diss, control$prop.args.diss, nw, class="d")
  model.diss <- ergm.getmodel(dissolution, nw, stergm.order=stergm.order, MHp=MHproposal.diss)
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
                  verbose)
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
  v$etamap <- model.form$etamap
  options(warn=current.warn)
  v
}
