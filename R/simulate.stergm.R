#========================================================================
# This file contains the following 2 functions for simulating stergms
#           <simulate.stergm>
#           <simulate.formula.stergm>
#========================================================================


########################################################################
# Each of the <simulate.X> functions collects a given number of networks
# drawn from the given distribution on the set of all networks; these
# may be returned as only the vector/matrix of sufficient statistics or
# as the networks and their statistics
#
# --PARAMETERS--
#   object       : either a stergm or a formation formula of the form
#                  'nw ~ term(s)'
#   dissolution  : for a formula 'object', this is the corresponding
#                  dissolution formula
#   nsim         : the number of networks to draw; default=1
#   seed         : an integer at which to set the random generator;
#                  default=NULL
#   theta.form   : the initial theta formation coefficients;
#                  default='object'$coef.form for stergm objects and
#                  default= a vector of 0's for formula objects
#   theta.diss   : the initial theta dissolution coefficients;
#                  default='object'$coef.diss for stergm objects and
#                  default= a vector of 0's for formula objects
#   time.burnin  : the number of MCMC steps to disregard before any MCMC
#                  sampling is done; default=0
#   time.interval: the number of MCMC steps between sampled networks;
#                  default=1
#   constraints  : a one-sided formula specifying the constraints on the
#                  support of the distribution of networks being simulated;
#                  default='object'$constraints for stergms, "~." for formulas
#   control      : a list of control parameters for algorithm tuning;
#                  default=<control.simulate.stergm>
#   toggles      : whether 'changed', the toggle matrix of timestamps and
#                  toggles, should be included in the return list (T or F);
#                  'toggles' will be switched to FALSE if either of
#                  'time.burnin' or 'time.interval' do not have their default
#                  values; default=TRUE
#   verbose      : whether to print out information on the status of
#                  the simulations; default=FALSE
#
# --RETURNED--
#   only 
#      nw:  the final network from the simulation routine
#   if 'control$final'=TRUE (the default is FALSE)
#   otherwise
#     outlist: a network.series object as a list containing:
#        formation  : the formation formula
#        dissolution: the dissolution formula
#        coef.form  : the passed in or defaulted 'coef.form'
#        coef.diss  : the passed in or defaulted 'coef.diss'
#        networks   : the list of simulated networks
#        constraints: the constraints formula
#        stats.form : the matrix of sampled statistics for 'model.form'
#        stats.diss : the matrix of sampled statistics for 'model.form'
#        changed    : a toggle matrix, where the first column is
#                     the timestamp of the toggle and the 2nd and 3rd
#                     columns are the head & tail of the toggle; this
#                     is only returned if the input param 'toggles'
#                     ends up being TRUE (see above) 
#                    'changed' will also have 2 attributes:
#            start  : 1
#            end    : the number of simulations
#        maxchanges : the size of "MCMC Dyn workspace"
#
###############################################################################

simulate.stergm<-function(object,
                          nsim=1,
                          coef.form=object$formation.fit$coef,coef.diss=object$dissolution.fit$coef,
                          time.burnin=0, time.interval=1,
                          control=control.simulate.stergm(),
                          statsonly=FALSE,
                          toggles=!statsonly,
                          verbose=FALSE, ...){
  simulate.formula.stergm(object$formation,dissolution=object$dissolution,nsim=nsim,coef.form=coef.form, coef.diss=coef.diss,  time.burnin=time.burnin, time.interval=time.interval,control=control,toggles=toggles,verbose=verbose,...)
}



# Note that we are overriding simulate.network here, since the first argument is a network.
simulate.network<-simulate.formula.stergm <- function(object, nsim=1,
                                                      formation, dissolution,
                                                      coef.form,coef.diss,
                                                      time.burnin=0, time.interval=1,
                                                      control=control.simulate.stergm(),
                                                      toggles = TRUE,
                                                      verbose=FALSE) {
  
  if(!is.null(control$seed)) set.seed(as.integer(control$seed))
  # Toggles is a "main call" parameter, since it affects what to
  # compute rather than just how to compute it, but it's convenient to
  # have it as a part of the control data structure.
  control$toggles <- toggles
  
  nw <- as.network(object)
  if(!is.network(nw)){
    stop("A network object on the LHS of the formation must be given")
  }

  formation<-ergm.update.formula(formation,nw~.)
  dissolution<-ergm.update.formula(dissolution,nw~.)
  
  model.form <- ergm.getmodel(formation, nw)
  if(!missing(coef.form) && coef.length.model(model.form)!=length(coef.form)) stop("coef.form has ", length(coef.form), " elements, while the model requires ",coef.length.model(model.form)," parameters.")

  model.diss <- ergm.getmodel(dissolution, nw)
  if(!missing(coef.diss) && coef.length.model(model.diss)!=length(coef.diss)) stop("coef.diss has ", length(coef.diss), " elements, while the model requires ",coef.length.model(model.diss)," parameters.")


  
  verbose <- match(verbose,
                c("FALSE","TRUE", "very"), nomatch=1)-1
  if(missing(coef.form)) {
    coef.form <- rep(0,length(model.form$coef.names))
    warning("No parameter values given, using Bernouli network.\n\t")
  }

  if(missing(coef.diss)) {
    coef.diss <- rep(0,length(model.diss$coef.names))
    warning("No parameter values given, using Bernoulli dissolution.\nThis means that every time step, half the ties get dissolved!\n\t")
  }

  if((time.burnin!=0 || time.interval!=1) && control$toggles){
    warning("Burnin is present or interval isn't 1. Toggle list will not be returned.")
    control$toggles<-FALSE
  }
    
  MHproposal.form <- MHproposal(~.,control$form$MCMC.prop.args,nw,
                                weights=control$form$MCMC.prop.weights,class="f")
  MHproposal.diss <- MHproposal(~.,control$diss$MCMC.prop.args,nw,
                                weights=control$diss$MCMC.prop.weights,class="d")

  eta.form <- ergm.eta(coef.form, model.form$etamap)
  eta.diss <- ergm.eta(coef.diss, model.diss$etamap)

  control$time.burnin <- time.burnin
  control$time.interval <- time.interval
  control$time.samplesize <- nsim
  
  z <- stergm.getMCMCsample(nw, model.form, model.diss,
                            MHproposal.form, MHproposal.diss,
                            eta.form, eta.diss, control, verbose)
  

  # FIXME: Standardize output format once the conversion function is available.
  if(toggles){
    changed<-z$changed
    attr(changed,"start")<-time.burnin+1
    attr(changed,"end")<-(nsim-1)*time.interval+time.burnin+1
    out.list <- list(formation = formation,
                     dissolution = dissolution,
                     networks = nw,
                     changed=changed, 
                     maxchanges=z$maxchanges,
                     stats.form = z$statsmatrix.form,stats.diss = z$statsmatrix.diss,
                     coef.form=coef.form,coef.diss=coef.diss)
    class(out.list) <- "network.series"
  }else{out.list<-list(stats.form = z$statsmatrix.form,stats.diss = z$statsmatrix.diss)}
  return(out.list)
}
