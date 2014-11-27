#  File R/simulate.ergm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
#========================================================================
# This file contains the following 2 functions for simulating ergms
#           <simulate.ergm>
#           <simulate.formula.ergm>
#========================================================================


########################################################################
# Each of the <simulate.X> functions collects a given number of networks
# drawn from the given distribution on the set of all networks; these
# may be returned as only the vector/matrix of sufficient statistics or
# as the networks and their statistics
#
# --PARAMETERS--
#   object     : either aern ergm or a formula of the form 'nw ~ term(s)'
#   nsim       : the number of networks to draw; default=1
#   basis      : optionally, a network to start the MCMC algorithm from;
#                if provided, this overrides the network given in
#                'object's formula; default=NULL
#   seed       : an integer at which to set the random generator;
#                default=NULL
#   coef     : the set of parameters from which the sample is to be
#                drawn; default='object$coef' if 'object' is an ergm or
#                0 if a formula
#   burnin     : the number of proposals to disregard before any MCMC
#                sampling is done; default=1000
#   interval   : the number of proposals between sampled networks;
#                default=1000
#   statsonly  : whether to return only the network statistics;
#                default=FALSE
#   sequential : whether subsequent draws should use the prior draw
#                as the starting network in the MCMC algorithm (T or F);
#                if FALSE, the initial network is always used as the
#                starting network; default=TRUE
#   constraints: a one-sided formula specifying the constraints on the
#                support of the distribution of networks being simulated;
#                default=NULL
#   control    : a list of control parameters for algorithm tuning, as
#                returned by <control.simulate.ergm> or
#                <control.simulate.formula>; default=<control.simulate.X>
#   verbose    : whether to print out information on the status of
#                the simulations; default=FALSE
#
# --RETURNED--
#   if 'statsonly'=TRUE  -- the vector of summary statistics for the
#      'nsim'=1             drawn network
#   if 'statsonly'=TRUE  -- the matrix of summary statistics for each
#                           drawn network; each row corresponds to a network
#   if 'statsonly'=FALSE -- the drawn network
#      'nsim'=1
#   if 'statsonly'=FALSE -- a list with the following components:
#      'nsim'>1              formula : 'object'
#                            networks: the list of drawn networks
#                            stats   : the matrix of summary stats
#                            coef    : 'init'
#
###############################################################################

simulate.ergm <- function(object, nsim=1, seed=NULL, 
                          coef=object$coef,
                          response=object$response,
                          reference=object$reference,
                          constraints=object$constraints,
                          monitor=NULL,
                          statsonly=FALSE,
                          esteq=FALSE,
                          sequential=TRUE,
                          control=control.simulate.ergm(),
                          verbose=FALSE, ...) {
  check.control.class(c("simulate.ergm","simulate.formula"))
  control.transfer <- c("MCMC.burnin", "MCMC.interval", "MCMC.prop.weights", "MCMC.prop.args", "MCMC.packagenames", "MCMC.init.maxedges")
  for(arg in control.transfer)
    if(is.null(control[[arg]]))
      control[arg] <- list(object$control[[arg]])

  control <- set.control.class("control.simulate.formula")
  
  simulate.formula(object$formula, nsim=nsim, coef=coef, response=response, reference=reference,
                   statsonly=statsonly,
                   esteq=esteq,
                   sequential=sequential, constraints=constraints,
                   monitor=monitor,
                   control=control, verbose=verbose, seed=seed, ...)
}


simulate.formula <- function(object, nsim=1, seed=NULL,
                               coef, response=NULL, reference=~Bernoulli,
                               constraints=~.,
                               monitor=NULL,
                               basis=NULL,
                               statsonly=FALSE,
                               esteq=FALSE,
                               sequential=TRUE,
                               control=control.simulate.formula(),
                               verbose=FALSE, ...) {
  check.control.class(myname="ERGM simulate.formula")
  # Backwards-compatibility code:
  if("theta0" %in% names(list(...))){
    warning("Passing the parameter vector as theta0= is deprecated. Use coef= instead.")
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
  form <- ergm.update.formula(object, basis ~ ., from.new="basis")

  if(!is.null(monitor)){
    # Construct a model to get the number of parameters monitor requires.
    monitor <- ergm.update.formula(monitor, nw~., from.new="nw")
    monitor.m <- ergm.getmodel(monitor, basis, response=response)
    monitored.length <- coef.length.model(monitor.m)
    
    monitor <- term.list.formula(monitor[[3]])
    form<-append.rhs.formula(form, monitor)
  }else{
    monitored.length <- 0
  }

  # Prepare inputs to ergm.getMCMCsample
  m <- ergm.getmodel(form, basis, response=response, role="static")
  # Just in case the user did not give a coef value, set it to zero.
  # (probably we could just return an error in this case!)
  if(missing(coef)) {
    coef <- c(rep(0, coef.length.model(m)))
    warning("No parameter values given, using Bernouli network\n\t")
  }

  coef <- c(coef, rep(0, monitored.length))
  
  if(coef.length.model(m)!=length(coef)) stop("coef has ", length(coef) - monitored.length, " elements, while the model requires ",coef.length.model(m) - monitored.length," parameters.")

  MHproposal <- MHproposal(constraints,arguments=control$MCMC.prop.args,
                           nw=nw, weights=control$MCMC.prop.weights, class="c",reference=reference,response=response)  

  if (any(is.nan(coef) | is.na(coef)))
    stop("Illegal value of coef passed to simulate.formula")
  
  # Create eta0 from coef
  eta0 <- ergm.eta(coef, m$etamap)
    
  # Create vector of current statistics
  curstats<-summary(form,response=response)
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
    # In this case, we can make one, parallelized run of
    # ergm.getMCMCsample.
    control$MCMC.samplesize <- nsim
    z <- ergm.getMCMCsample(nw, m, MHproposal, eta0, control, verbose=verbose, response=response)
    
    # Post-processing:  Add term names to columns and shift each row by
    # observed statistics.
    colnames(z$statsmatrix) <- m$coef.names
    out.mat <- sweep(z$statsmatrix[seq_len(nsim),,drop=FALSE], 2, curstats, "+")
  }else{
    # Create objects to store output
    if (!statsonly) { 
      nw.list <- list()
    }
    out.mat <- matrix(nrow=0, ncol=length(curstats), 
                      dimnames = list(NULL, m$coef.names)) 
    
    # Call ergm.getMCMCsample once for each network desired.  This is much slower
    # than when sequential==TRUE and statsonly==TRUE, but here we have a 
    # more complicated situation:  Either we want a network for each
    # MCMC iteration (statsonly=FALSE) or we want to restart each chain
    # at the original network (sequential=FALSE).
    if(control$parallel) curstats <- matrix(curstats, nrow=control$parallel, ncol=length(curstats), byrow=TRUE)
    
    for(i in 1:ceiling(nsim/max(control$parallel,1))){
      
      control$MCMC.samplesize <- if(control$parallel==0) 1 else control$parallel
      control$MCMC.burnin <- if(i==1 || sequential==FALSE) control$MCMC.burnin else control$MCMC.interval
      z <- ergm.getMCMCsample(nw, m, MHproposal, eta0, control, verbose=verbose, response=response)
      
      out.mat <- rbind(out.mat, curstats + z$statsmatrix)
      
      if(!statsonly) # then store the returned network:
        if(control$parallel==0) nw.list[[length(nw.list)+1]] <- z$newnetwork else nw.list <- c(nw.list, z$newnetworks)
      
      if(sequential){ # then update the network state:
        nw <- if(control$parallel==0) z$newnetwork else z$newnetworks
        curstats <- curstats + z$statsmatrix
      }

      if(verbose){cat(sprintf("Finished simulation %d of %d.\n",i, nsim))}
    }
  }

  out.mat <- out.mat[seq_len(nsim),,drop=FALSE]

  if(esteq) out.mat <- .ergm.esteq(coef, m, out.mat)
  
  if (statsonly)
    return(out.mat)
  
  # If we get here, statsonly==FALSE.
  if (nsim==1) {
    return(nw.list[[1]])
  } else {
    nw.list <- nw.list[seq_len(nsim)]
    attributes(nw.list) <- list(formula=object, stats=out.mat, coef=coef,
                                control=control,
                                constraints=constraints, reference=reference,
                                monitor=monitor, response=response)

    class(nw.list) <- "network.list"
    return(nw.list)
  }
}


