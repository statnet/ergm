#  File ergm/R/simulate.stergm.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
########################################################################
# Each of the <simulate.X> functions collects a given number of networks
# drawn from the given distribution on the set of all networks; these
# may be returned as only the vector/matrix of sufficient statistics or
# as the networks and their statistics
###############################################################################

simulate.stergm<-function(object, nsim=1, seed=NULL,
                          coef.form=object$formation.fit$coef,coef.diss=object$dissolution.fit$coef,
                          monitor = object$targets,
                          time.slices, time.burnin=0, time.interval=1,
                          control=control.simulate.stergm(),
                          statsonly=time.burnin>0||time.interval>1,
                          stats.form = FALSE,
                          stats.diss = FALSE,
                          verbose=FALSE, ...){
  simulate.network(object$network,formation=object$formation,dissolution=object$dissolution,nsim=nsim,coef.form=coef.form, coef.diss=coef.diss, monitor=monitor, time.slices=time.slices, time.burnin=time.burnin, time.interval=time.interval,control=control, statsonly=statsonly, stats.form = stats.form, stats.diss = stats.diss, verbose=verbose,...)
}



# Note that we are overriding simulate.network here, since the first argument is a network.
simulate.network <- function(object, nsim=1, seed=NULL,
                             formation, dissolution,
                             coef.form,coef.diss,
                             monitor = NULL,
                             time.slices, time.burnin=0, time.interval=1,
                             control=control.simulate.stergm(),
                             statsonly=time.burnin>0||time.interval>1,
                             stats.form = FALSE,
                             stats.diss = FALSE,
                             verbose=FALSE, ...) {
  if(length(list(...))) stop("Unknown arguments: ",names(list(...)))
  
  if(!is.null(seed)) set.seed(as.integer(seed))
  # Toggles is a "main call" parameter, since it affects what to
  # compute rather than just how to compute it, but it's convenient to
  # have it as a part of the control data structure.
  control$toggles <- !statsonly

  nw <- as.network(object)
  if(!is.network(nw)){
    stop("A network object must be given")
  }

  formation<-ergm.update.formula(formation,nw~.)
  dissolution<-ergm.update.formula(dissolution,nw~.)

  unset.offset.call <- function(object){
    if(inherits(object,"call") && object[[1]]=="offset")
      object[[2]]
    else
      object
  }
  
  if(is.character(monitor)){
    monitor <- switch(monitor,
                      formation = formation,
                      dissolution = dissolution,
                      all = append.rhs.formula(~nw, unique(lapply(c(term.list.formula(formation[[3]]),term.list.formula(dissolution[[3]])), unset.offset.call)))
                      )
  }
  
  if(!is.null(monitor)) monitor<-ergm.update.formula(monitor,nw~.)
  
  model.form <- ergm.getmodel(formation, nw, role="formation")
  if(!missing(coef.form) && coef.length.model(model.form)!=length(coef.form)) stop("coef.form has ", length(coef.form), " elements, while the model requires ",coef.length.model(model.form)," parameters.")

  model.diss <- ergm.getmodel(dissolution, nw, role="dissolution")
  if(!missing(coef.diss) && coef.length.model(model.diss)!=length(coef.diss)) stop("coef.diss has ", length(coef.diss), " elements, while the model requires ",coef.length.model(model.diss)," parameters.")

  model.mon <- if(!is.null(monitor)) ergm.getmodel(monitor, nw, role="target") else NULL
  
  verbose <- match(verbose,
                c("FALSE","TRUE", "very"), nomatch=1)-1
  if(missing(coef.form)) {
    coef.form <- rep(0,length(model.form$coef.names))
    warning("No parameter values given, using Bernouli formation.\nThis means that every time step, half the non-tie dyads will gain a tie!")
  }

  if(missing(coef.diss)) {
    coef.diss <- rep(0,length(model.diss$coef.names))
    warning("No parameter values given, using Bernoulli dissolution.\nThis means that every time step, half the ties get dissolved!\n")
  }

  if((time.burnin!=0 || time.interval!=1) && control$toggles){
    warning("Burnin is present or interval isn't 1. Toggle list will not be returned.")
    control$toggles<-FALSE
  }
    
  MHproposal.form <- MHproposal(~.,control$MCMC.prop.args.form,nw,
                                weights=control$MCMC.prop.weights.form,class="f")
  MHproposal.diss <- MHproposal(~.,control$MCMC.prop.args.diss,nw,
                                weights=control$MCMC.prop.weights.diss,class="d")

  eta.form <- ergm.eta(coef.form, model.form$etamap)
  eta.diss <- ergm.eta(coef.diss, model.diss$etamap)
  
  control$time.burnin <- time.burnin
  control$time.interval <- time.interval
  control$time.samplesize <- time.slices
  control$collect.form <- stats.form
  control$collect.diss <- stats.diss

  out <- replicate(nsim, {
    if(is.null(nw %n% "lasttoggle")) nw %n% "lasttoggle" <- rep(0, network.dyadcount(nw))
    if(is.null(nw %n% "time")) nw %n% "time" <- 0
    
    z <- stergm.getMCMCsample(nw, model.form, model.diss, model.mon,
                              MHproposal.form, MHproposal.diss,
                              eta.form, eta.diss, control, verbose)
    
    stats.form <- if(control$collect.form) mcmc(sweep(z$statsmatrix.form,2,summary(formation),"+"),start=time.burnin+1,thin=time.interval)
    stats.diss <- if(control$collect.diss) mcmc(sweep(z$statsmatrix.diss,2,summary(dissolution),"+"),start=time.burnin+1,thin=time.interval)
    stats <- if(!is.null(model.mon)) mcmc(sweep(z$statsmatrix.mon,2,summary(monitor),"+"),start=time.burnin+1,thin=time.interval)
    
    if(control$toggles){
      library(networkDynamic)
      nwd <- as.networkDynamic(nw, toggles = z$changed, start = nw%n%"time" + 0, end = nw%n%"time" + time.slices)
      nwd<-delete.network.attribute(nwd, "time")
      nwd<-delete.network.attribute(nwd, "lasttoggle")
      attributes(nwd) <- c(attributes(nwd), # Don't clobber existing attributes!
                           list(formation = formation,
                                dissolution = dissolution,
                                stats.form = stats.form,
                                stats.diss = stats.diss,
                                stats = stats,
                                coef.form=coef.form,
                                coef.diss=coef.diss,
                                start = nw%n%"time" + 0,
                                end = nw%n%"time" + time.slices,
                                toggles = z$changed))
      nwd
    }else
    list(stats.form = stats.form,stats.diss = stats.diss, stats = stats)
  },
                   simplify = FALSE)
  if(nsim==1){
    out<-out[[1]]
    if(!control$toggles){
      for(name in names(out)) # Strip the unreturned stats matrices.
        if(is.null(out[[name]]))
          out[[name]] <- NULL
      if(length(out))
        out <- out[[1]] # If there is only one, just return it.
    }
    out
  }else{
    if(control$toggles){
      # If we've returned a list of networkDynamics, then it's a network list.
      class(out) <- "network.list"
      out
    }else{
      # Otherwise, we can combine the simulation into mcmc.lists.
      outl <- list()
      for(name in names(out[[1]])){
        if(!is.null(out[[1]][[name]]))
          outl[[name]] <- do.call(mcmc.list, lapply(out, "[[", name))
      }
      if(length(outl))
         outl <- outl[[1]] # If there is only one, just return it.
      outl
    }
  }
}
