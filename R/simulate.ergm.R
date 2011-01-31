#========================================================================
# This file contains the following 2 functions for simulating ergms
#           <simulate.ergm>
#           <simulate.formula>
#========================================================================


########################################################################
# Each of the <simulate.X> functions collects a given number of networks
# drawn from the given distribution on the set of all networks; these
# may be returned as only the vector/matrix of sufficient statistics or
# as the networks and their statistics
#
# --PARAMETERS--
#   object     : either an ergm or a formula of the form 'nw ~ term(s)'
#   nsim       : the number of networks to draw; default=1
#   basis      : optionally, a network to start the MCMC algorithm from;
#                if provided, this overrides the network given in
#                'object's formula; default=NULL
#   seed       : an integer at which to set the random generator;
#                default=NULL
#   theta0     : the set of parameters from which the sample is to be
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
#                            coef    : 'theta0'
#
###############################################################################

simulate.ergm <- function(object, nsim=1, seed=NULL, theta0=object$coef,
                          burnin=1000, interval=1000,
                          statsonly=FALSE,
                          sequential=TRUE,
                          constraints=NULL,
                          control=control.simulate.ergm(),
                          verbose=FALSE, ...) {
  if(is.null(burnin)){burnin <- object$burnin}
  if(is.null(interval)){interval <- object$interval}
  if(is.null(constraints)){constraints <- object$constraints}
  simulate.formula(object$formula, nsim=nsim, seed=seed, theta0=theta0,
                   burnin=burnin, interval=interval, statsonly=statsonly,
                   sequential=sequential, constraints=constraints,
                   control=control, verbose=verbose, ...)
}



simulate.formula <- function(object, nsim=1, seed=NULL, theta0,
                             burnin=1000, interval=1000,
                             basis=NULL,
                             statsonly=FALSE,
                             sequential=TRUE,
                             constraints=~.,
                             control=control.simulate.formula(),
                             verbose=FALSE, ...) {
  if(!is.null(seed)) {set.seed(as.integer(seed))}

  # If basis is not null, replace network in formula by basis.
  # In either case, let nw be network object from formula.
  if(is.null(nw <- basis)) {
    nw <- ergm.getnetwork(object)
  }
  
  # Do some error-checking on the nw object
  if(class(nw) =="network.series"){
    nw <- nw$networks[[1]]
  }
  nw <- as.network(nw)
  if(!is.network(nw)){
    stop("A network object on the LHS of the formula or via",
         " the 'basis' argument must be given")
  }

  # New formula (no longer use 'object'):
  form <- ergm.update.formula(object, nw ~ .)
  
  # Prepare inputs to ergm.getMCMCsample
  m <- ergm.getmodel(form, nw, drop=FALSE)
  Clist <- ergm.Cprepare(nw, m)
  MHproposal <- MHproposal(constraints,arguments=control$prop.args,
                           nw=nw, model=m, weights=control$prop.weights, class="c")  

  # Just in case the user did not give a theta0 value, set it to zero.
  # (probably we could just return an error in this case!)
  if(missing(theta0)) {
    theta0 <- rep(0,Clist$nstats)
    warning("No parameter values given, using Bernouli network\n\t")
  }
  if (any(is.infinite(theta0))){
   theta0[is.infinite(theta0)] <- sign(theta0[is.infinite(theta0)])*10000 
  }
  if (any(is.nan(theta0) | is.na(theta0)))
    stop("Illegal value of theta0 passed to simulate.formula")
    
  # Create vector of current statistics
  curstats<-summary(form)
  names(curstats) <- m$coef.names

  # prepare MCMCparams object
  MCMCparams <- list(samplesize=1,
                     maxedges = 1+max(control$maxedges, Clist$nedges),
#                     stats=curstats, # deprecated as of version 2.2-3
                     burnin=burnin,
                     interval=interval,
                     parallel=control$parallel,
                     packagenames=control$packagenames,
                     Clist.miss=ergm.design(nw, m, verbose=verbose))
  
  # Explain how many iterations and steps will ensue if verbose==TRUE
  if (verbose) {
    cat (paste ("Starting MCMC iterations to generate ", nsim,
                " network", ifelse(nsim>1,"s",""), sep=""))
  }
  
  #########################
  ## Main part of function:
  if(sequential && statsonly){ 
    # Call ergm.getMCMCsample only one time, using the C function to generate the whole
    # matrix of network statistics.
    MCMCparams$samplesize <- nsim
    MCMCparams$nmatrixentries <- nsim * length(curstats)
    z <- ergm.getMCMCsample(Clist, MHproposal, theta0, MCMCparams, verbose=verbose)
    
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
  MCMCparams$nmatrixentries <- length(curstats)
  for(i in 1:nsim){
    MCMCparams$burnin <- ifelse(i==1, burnin, interval)
    z <- ergm.getMCMCsample(Clist, MHproposal, theta0, MCMCparams, verbose)

    # Create a network object if statsonly==FALSE
    if (!statsonly) {
      nw.list[[i]] <- network.update(nw, z$newedgelist, matrix.type="edgelist",
                                     output=control$network.output)
    }
    out.mat[i,] <- curstats + z$statsmatrix
    if (sequential){ 
      # If we get here, statsonly must be FALSE
      nw <- as.network.uncompressed(nw.list[[i]])
      Clist <- ergm.Cprepare(nw, m)
      curstats <- curstats + z$statsmatrix
    }
    if(verbose){cat(sprintf("Finished simulation %d of %d.\n",i, nsim))}
  }
  
  if (statsonly)
    return(out.mat[1:nsim,]) # If nsim==1, this will return a vector, not a matrix
  
  # If we get here, statsonly==FALSE.
  if (nsim==1) {
    return(nw.list[[1]])
  } else {  
    out.list <- list(formula = object, networks = nw.list, 
                     stats = out.mat, coef=theta0)
    class(out.list) <- "network.series"
    return(out.list)
  }
}


