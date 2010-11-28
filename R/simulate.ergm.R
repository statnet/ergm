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
  simulate.formula(object$formula, nsim=nsim, seed=seed, theta0=theta0, response=object$response, reference=if(is.null(object$reference)) "Bernoulli" else object$reference,
                   burnin=burnin, interval=interval, statsonly=statsonly,
                   sequential=sequential, constraints=constraints,
                   control=control, verbose=verbose, ...)
}

# In CRAN version, the following function is called simulate.formula
# Here, there is a good reason to call it simulate.formula.ergm:
# see simulate.formula.R
simulate.formula.ergm <- function(object, nsim=1, seed=NULL, theta0, response=NULL, reference="Bernoulli",
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
  formula <- ergm.update.formula(object, nw ~ .)
  
  # Prepare inputs to ergm.getMCMCsample
  m <- ergm.getmodel(formula, nw, drop=FALSE, response=response)
  Clist <- ergm.Cprepare(nw, m, response=response)
  MHproposal <- MHproposal(constraints,arguments=control$prop.args,
                           nw=nw, model=m, weights=control$prop.weights, class="c",reference=reference)  

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
  curstats<-summary(formula,response=response)
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
    cat(paste("Starting ",nsim," MCMC iteration", ifelse(nsim>1,"s",""),
              " of ", burnin+interval*(MCMCparams$samplesize-1), 
              " steps", ifelse(nsim>1, " each", ""), ".\n", sep=""))
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
  if(sequential){
   for(i in 1:nsim){
    MCMCparams$burnin <- ifelse(i==1 || !sequential, burnin, interval)
    z <- ergm.getMCMCsample(Clist, MHproposal, theta0, MCMCparams, verbose)

    # Create a network object if statsonly==FALSE
    if (!statsonly) {
      nw.list[[i]] <- newnw.extract(nw, z, output=control$network.output, response=response)
    }
    out.mat[i,] <- curstats + z$statsmatrix
    if (sequential){ 
      # If we get here, statsonly must be FALSE
      nw <- as.network.uncompressed(nw.list[[i]])
      Clist <- ergm.Cprepare(nw, m, response=response)
      curstats <- curstats + z$statsmatrix
    }
    if(verbose){cat(sprintf("Finished simulation %d of %d.\n",i, nsim))}
   }
  }else{
    MCMCparams.parallel <- MCMCparams
    MCMCparams.parallel$samplesize <- 1
    MCMCparams.parallel$nmatrixentries <- length(curstats)
    require(snow)
#
# Start PVM if necessary
#
    if(getClusterOption("type")=="PVM"){
     if(verbose){cat("Engaging warp drive using PVM ...\n")}
     require(rpvm)
     PVM.running <- try(.PVM.config(), silent=TRUE)
     if(inherits(PVM.running,"try-error")){
      hostfile <- paste(Sys.getenv("HOME"),"/.xpvm_hosts",sep="")
      .PVM.start.pvmd(hostfile)
      cat("no problem... PVM started by ergm...\n")
     }
    }else{
     if(verbose){cat("Engaging warp drive using MPI ...\n")}
    }
#
#   Start Cluster
#
    cl<-makeCluster(MCMCparams$parallel)
    clusterSetupRNG(cl)
    if("ergm" %in% MCMCparams$packagenames){
     clusterEvalQ(cl,library(ergm))
    }
#
#   Run the jobs with rpvm or Rmpi
#
    flush.console()
    outlist <- clusterCall(cl,ergm.getMCMCsample,
     Clist,MHproposal,theta0,MCMCparams.parallel,verbose)
#
#   Process the results
#
    for(i in (1:MCMCparams$parallel)){
     z <- outlist[[i]]
     nedges <- z$newnwheads[1]
     # Create a network object if statsonly==FALSE
     if (!statsonly) {
      nw.list[[i]] <- newnw.extract(nw, z, output=control$network.output, response=response)
     }
     nw <- as.network.uncompressed(nw.list[[i]])
     out.mat[i,] <- summary(formula)
    }
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


