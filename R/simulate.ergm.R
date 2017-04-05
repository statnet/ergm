#  File ergm/R/simulate.ergm.R
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

simulate.ergm <- function(object, nsim=1, seed=NULL, 
                          coef=object$coef,
                          constraints=object$constraints,
                          monitor=NULL,
                          statsonly=FALSE,
                          sequential=TRUE,
                          control=control.simulate.ergm(),
                          verbose=FALSE, ...) {
  control.transfer <- c("MCMC.burnin", "MCMC.interval", "MCMC.prop.weights", "MCMC.prop.args", "MCMC.packagenames", "MCMC.init.maxedges")
  for(arg in control.transfer)
    if(is.null(control[[arg]]))
      control[arg] <- list(object$control[[arg]])

  simulate.formula(object$formula, nsim=nsim, coef=coef,
                   statsonly=statsonly,
                   sequential=sequential, constraints=constraints,
                   monitor=monitor,
                   control=control, verbose=verbose, seed=seed, ...)
}


simulate.formula <- function(object, nsim=1, seed=NULL,
                               coef, 
                               constraints=~.,
                               monitor=NULL,
                               basis=NULL,
                               statsonly=FALSE,
                               sequential=TRUE,
                               control=control.simulate.formula(),
                               verbose=FALSE, ...) {
  # Backwards-compatibility code:
  if("theta0" %in% names(list(...))){
    warning("Passing the parameter vector as theta0= is depcrecated. Use coef= instead.")
    coef<-list(...)$theta0
  }
  control <- control.simulate.ergm.toplevel(control,...)
  
  if(!is.null(seed)) {set.seed(as.integer(seed))}
  
  # define nw as either the basis argument or (if NULL) the LHS of the formula
  if (is.null(nw <- basis)) {
    nw <- ergm.getnetwork(object)    
  }
  
  # Do some error-checking on the nw object
  nw <- as.network(nw)
  if(!is.network(nw)){
    stop("A network object on the LHS of the formula or via",
         " the 'basis' argument must be given")
  }
  if(is.null(basis)) {
    basis <- nw
  }

  # New formula (no longer use 'object'):
  form <- ergm.update.formula(object, basis ~ .)

  if(!is.null(monitor)){
    # Construct a model to get the number of parameters monitor requires.
    monitor <- ergm.update.formula(monitor, nw~.)
    monitor.m <- ergm.getmodel(monitor, basis)
    monitored.length <- coef.length.model(monitor.m)
    
    monitor <- term.list.formula(monitor[[3]])
    form<-append.rhs.formula(form, monitor)
  }else{
    monitored.length <- 0
  }

  # Prepare inputs to ergm.getMCMCsample
  m <- ergm.getmodel(form, basis, role="static")
  # Just in case the user did not give a coef value, set it to zero.
  # (probably we could just return an error in this case!)
  if(missing(coef)) {
    coef <- c(rep(0, coef.length.model(m)))
    warning("No parameter values given, using Bernouli network\n\t")
  }

  coef <- c(coef, rep(0, monitored.length))
  
  if(coef.length.model(m)!=length(coef)) stop("coef has ", length(coef) - monitored.length, " elements, while the model requires ",coef.length.model(m) - monitored.length," parameters.")

  MHproposal <- MHproposal(constraints,arguments=control$MCMC.prop.args,
                           nw=nw, weights=control$MCMC.prop.weights, class="c")  

  if (any(is.infinite(coef))){
   coef[is.infinite(coef)] <- sign(coef[is.infinite(coef)])*10000 
  }
  if (any(is.nan(coef) | is.na(coef)))
    stop("Illegal value of coef passed to simulate.formula")
  
  # Create eta0 from coef
  eta0 <- ergm.eta(coef, m$etamap)
    
  # Create vector of current statistics
  curstats<-summary(form)
  names(curstats) <- m$coef.names

  # prepare control object
  control$MCMC.init.maxedges <- 1+max(control$MCMC.init.maxedges, network.edgecount(nw))
  
  # Explain how many iterations and steps will ensue if verbose==TRUE
  if (verbose) {
    cat (paste ("Starting MCMC iterations to generate ", nsim,
                " network", ifelse(nsim>1,"s\n","\n"), sep=""))
  }
  
  #########################
  ## Main part of function:
  if(sequential && statsonly){ 
    # Call ergm.getMCMCsample only one time, using the C function to generate the whole
    # matrix of network statistics.
    control$MCMC.samplesize <- nsim
    z <- ergm.getMCMCsample(nw, m, MHproposal, eta0, control, verbose=verbose)
    
    # Post-processing:  Add term names to columns and shift each row by
    # observed statistics.
    colnames(z$statsmatrix) <- m$coef.names
    return(sweep(z$statsmatrix, 2, curstats, "+"))
  } 
  
  # If we get here, either sequential==FALSE or statsonly==FALSE.
  # Create objects to store output
  if (!statsonly) { 
    nw.list <- list()
  }
  out.mat <- matrix(nrow=nsim, ncol=length(curstats), 
                    dimnames = list(NULL, m$coef.names)) 
  
  # Call ergm.getMCMCsample once for each network desired.  This is much slower
  # than when sequential==TRUE and statsonly==TRUE, but here we have a 
  # more complicated situation:  Either we want a network for each
  # MCMC iteration (statsonly=FALSE) or we want to restart each chain
  # at the original network (sequential=FALSE).
  if (sequential) { # non-parallel method used here
    for(i in 1:nsim){
      control$MCMC.samplesize <- 1
      control$MCMC.burnin <- ifelse(i==1, control$MCMC.burnin, control$MCMC.interval)
      z <- ergm.getMCMCsample(nw, m, MHproposal, eta0, control, verbose=verbose)
      
      # Create a network object if statsonly==FALSE
      if (!statsonly) {
        nw <- nw.list[[i]] <- z$newnetwork
      }
      out.mat[i,] <- curstats + z$statsmatrix
      # If we get here, statsonly must be FALSE
      curstats <- curstats + z$statsmatrix
      if(verbose){cat(sprintf("Finished simulation %d of %d.\n",i, nsim))}
    }
  } else {
    # non-sequential
    control$MCMC.samplesize <- 1
    z <- ergm.getMCMCsample(nw, m, MHproposal, eta0, control, verbose=verbose)
    if (!statsonly) {
      nw.list <- z$newnetworks
    }

  }
  
  if (statsonly)
    return(out.mat[1:nsim,]) # If nsim==1, this will return a vector, not a matrix
  
  # If we get here, statsonly==FALSE.
  if (nsim==1) {
    return(nw.list[[1]])
  } else {
    attributes(nw.list) <- list(formula=object, stats=out.mat, coef=coef,
                                control=control,
                                constraints=constraints,
                                monitor=monitor)

    class(nw.list) <- "network.list"
    return(nw.list)
  }
}


