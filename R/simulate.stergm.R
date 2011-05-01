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
#    MH.burnin   : this is received as 'MH_interval' and determines the
#                  number of MH proposals used in each MCMC step;
#                  default='object'$MH.burnin for stergm objects;
#                  default=1000 for formula objects
#   constraints  : a one-sided formula specifying the constraints on the
#                  support of the distribution of networks being simulated;
#                  default='object'$constraints for stergms, "~." for formulas
#   stergm.order : the string describing the formation and dissolution order;
#                  default='object'$stergm.order for stergms, "DissAndForm"
#                  for formulas
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
                          nsim=1, seed=NULL,
                          theta.form=object$coef.form,theta.diss=object$coef.diss,
                        time.burnin=0, time.interval=1, MH.burnin=object$MH.burnin,
                        constraints=object$constraints,
                        stergm.order=object$stergm.order,
                        control=control.simulate.stergm(),
                        toggles=TRUE,
                        verbose=FALSE, ...){
  simulate.formula.stergm(object$formation,object$dissolution,nsim=nsim,seed=seed,theta.form=theta.form, theta.diss=theta.diss,  time.burnin=time.burnin, time.interval=time.interval,MH.burnin=MH.burnin,constraints=constraints,stergm.order=stergm.order,control=control,verbose=verbose,...)
}




simulate.formula.stergm <- function(object, nsim=1, seed=NULL, ..., dissolution,
                                theta.form,theta.diss,
                        time.burnin=0, time.interval=1, MH.burnin=1000,
                        constraints=~.,
                        stergm.order="DissAndForm",
                        control=control.simulate.stergm(),
                        toggles=TRUE,
                        verbose=FALSE) {

  formation <- object
  
  if(!is.null(seed)) set.seed(as.integer(seed))
  nw <- ergm.getnetwork(formation)
  if(class(nw) =="network.series"){
    nw <- nw$networks[[1]]
  }
  nw <- as.network(nw)
  if(!is.network(nw)){
    stop("A network object on the LHS of the formation must be given")
  }

  dissolution<-ergm.update.formula(dissolution,nw~.)
  
  model.form <- ergm.getmodel(formation, nw, drop=control$drop)
  model.diss <- ergm.getmodel(dissolution, nw, stergm.order=stergm.order,drop=control$drop)

  verbose <- match(verbose,
                c("FALSE","TRUE", "very"), nomatch=1)-1
  if(missing(theta.form)) {
    theta.form <- rep(0,length(model.form$coef.names))
    warning("No parameter values given, using Bernouli network.\n\t")
  }

  if(missing(theta.diss)) {
    theta.diss <- rep(0,length(model.diss$coef.names))
    warning("No parameter values given, using Bernoulli dissolution.\nThis means that every time step, half the ties get dissolved!\n\t")
  }

  if((time.burnin!=0 || time.interval!=1) && toggles){
    warning("Burnin is present or interval isn't 1. Toggle list will not be returned.")
    toggles<-FALSE
  }
  
  if(!is.null(seed)) set.seed(as.integer(seed))
    
  MHproposal.form <- MHproposal(constraints,control$prop.args.form,nw,
                                                    model.form,weights=control$prop.weights.form,class="f")
  MHproposal.diss <- MHproposal(constraints,control$prop.args.diss,nw,
                                                    model.diss,weights=control$prop.weights.diss,class="d")
  MCMCparams <- c(control,list(samplesize=nsim, time.interval=time.interval,
                               time.burnin=time.burnin,
                               MH.burnin=MH.burnin,
                               parallel=0,
                               meanstats.form=summary(model.form$formula),
                               meanstats.diss=summary(model.diss$formula),
                               toggles=toggles))
  
  z <- stergm.getMCMCsample(nw, model.form, model.diss,
                             MHproposal.form, MHproposal.diss,
                             theta.form, theta.diss, MCMCparams, verbose)

  if(control$final){
   nw <- z$newnetwork
   return(nw)
  }else{
    changed<-z$changed
    attr(changed,"start")<-time.burnin+1
    attr(changed,"end")<-(nsim-1)*time.interval+time.burnin+1
    out.list <- list(formation = formation,
                     dissolution = dissolution,
                     networks = nw,
                     constraints=constraints,
                     changed=changed, 
                     maxchanges=z$maxchanges,
                     stats.form = z$statsmatrix.form,stats.diss = z$statsmatrix.diss,
                     coef.form=theta.form,coef.diss=theta.diss)
    class(out.list) <- "network.series"
    return(out.list)
  }
}
