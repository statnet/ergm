#==============================================================================
# This file contains the 2 following functions for getting an MCMC sample
#      <ergm.getMCMCsample.parallel>
#      <ergm.mcmcslave>
#==============================================================================




#########################################################################################
# The <ergm.getMCMCsample.parallel> function samples networks using an MCMC algorithm via
# <MCMC_wrapper.C>. Unlike its <ergm.getMCMCsample> counterpart, this function is
# caple of running in multiple threads.  Note that the returned stats will be relative to
# the original network, i.e., the calling function must shift the statistics if required. 
# The calling function must also attach column names to the statistics matrix if required.
#
# --PARAMETERS--
#   nw        :  a network object
#   model     :  a model for the given 'nw' as returned by <ergm.getmodel>
#   MHproposal:  a list of the parameters needed for Metropolis-Hastings proposals and
#                the result of calling <MHproposal>
#   eta0      :  the initial eta coefficients
#   verbose   :  whether the C functions should be verbose; default=FALSE
#   MCMCparams:  list of MCMC tuning parameters; those recognized include
#       parallel    : the number of threads in which to run the sampling
#       packagenames: names of packages; this is only relevant if "ergm" is given
#       samplesize  : the number of networks to be sampled
#       interval    : the number of proposals to ignore between sampled networks
#       burnin      : the number of proposals to initially ignore for the burn-in
#                     period
#
# Note:  In reality, there should be many fewer arguments to this function,
# since most info should be passed via Clist (this is, after all, what Clist
# is for:  Holding all arguments required for the .C call).  In particular,
# the elements of MHproposal, MCMCparams, verbose should certainly
# be part of Clist.  But this is a project for another day!
#
# --RETURNED--
#   the sample as a list containing:
#     statsmatrix:  the stats matrix for the sampled networks, RELATIVE TO THE ORIGINAL
#                   NETWORK!
#     newnetwork :  the edgelist of the final sampled network
#     meanstats  :  NULL always (I think - 'meanstats' is taken from Clist, which
#                   is taken from a call to <ergm.Cprepare> which doesn't return
#                   a 'meanstats' component. Code to adjust 'meanstats' if null has
#                   been commented out)
#     nedges     :  the number of edges in the 'newnetwork' ??
#
#########################################################################################

ergm.getMCMCsample.parallel <- function(nw, model, MHproposal, eta0, MCMCparams, 
                                        verbose, response=NULL) {
  
  Clist <- ergm.Cprepare(nw, model, response=response)
  # maxedges <- max(5000, Clist$nedges)
  maxedges <- MCMCparams$maxedges/10

  # Run the first time, or as long as we have too many edges
  z <- NULL
  while(is.null(z) || z$status==1){
    maxedges <- 10*maxedges
    #
    #  Parallel running
    #
    if(MCMCparams$parallel==0){
      flush.console()
      z <- ergm.mcmcslave(Clist,MHproposal,eta0,MCMCparams,maxedges,verbose)
      if(z$status == 1){ # MCMC_TOO_MANY_EDGES
        warning("Sampled network has more edges than can be returned.")
        next
      }
      else if(z$status == 2){ # MCMC_MH_FAILED
        # MH proposal failed somewhere. Throw an error.
        error("Sampling failed due to a Metropolis-Hastings proposal failing.")
      }

      if(!is.null(z$burnin.failed) && z$burnin.failed) warning("Burn-in failed to converge after retries.")
      
      statsmatrix <- matrix(z$s, nrow=MCMCparams$samplesize,
                            ncol=Clist$nstats,
                            byrow = TRUE)
      newnetwork <- newnw.extract(nw,z,response=response)
    }else{
    MCMCparams.parallel <- MCMCparams
    MCMCparams.parallel$samplesize <- round(MCMCparams$samplesize / MCMCparams$parallel)

    cl <- ergm.getCluster(MCMCparams, verbose)
#
#   Run the jobs with rpvm or Rmpi
#
    flush.console()
    outlist <- clusterCall(cl,ergm.mcmcslave,
     Clist,MHproposal,eta0,MCMCparams.parallel,maxedges,verbose)
#
#   Process the results
#
    statsmatrix <- NULL
#   newedgelist <- matrix(0, ncol=2, nrow=0)
    for(i in (1:MCMCparams$parallel)){
     z <- outlist[[i]]
     if(z$status == 1){ # MCMC_TOO_MANY_EDGES
       warning("Sampled network has more edges than can be returned.")
       next
     }
     else if(z$status == 2){ # MCMC_MH_FAILED
       # MH proposal failed somewhere. Throw an error.
       error("Sampling failed due to a Metropolis-Hastings proposal failing.")
     }

     if(!is.null(z$burnin.failed) && z$burnin.failed) warning("Burn-in failed to converge after retries.")
     
     statsmatrix <- rbind(statsmatrix,
       matrix(z$s, nrow=MCMCparams.parallel$samplesize,
       ncol=Clist$nstats,
       byrow = TRUE))
#    if(z$newnw[1]>1){
#      newedgelist <- rbind(newedgelist,
     #                           matrix(z$newnw[2:z$newnw[1]], ncol=2, byrow=TRUE))
   }
    if(z$status == 1){ # MCMC_TOO_MANY_EDGES
      warning("Sampled network has more edges than can be returned.")
      next
    }
    newnetwork<-newnw.extract(nw,z,response=response)
    if(verbose){cat("parallel samplesize=",nrow(statsmatrix),"by",
	MCMCparams.parallel$samplesize,"\n")}

    ergm.stopCluster(cl)
  }
  }
  colnames(statsmatrix) <- model$coef.names
#
  statsmatrix[is.na(statsmatrix)] <- 0

##
## recenter statsmatrix by mean statistics if necessary
##
#   ms <- Clist$meanstats
#   if(!is.null(ms)) {
#     if (is.null(names(ms)) && length(ms) == length(model$coef.names))
#       names(ms) <- model$coef.names
##    obs <- summary(model$formula)
#     obs <- Clist$obs
##    print(paste("obs=",obs))
##    print(paste("statsmatrix=",apply(statsmatrix,2,mean)))
#     obs <- obs[match(colnames(statsmatrix), names(obs))]
#     ms  <-  ms[match(names(obs), names(ms))]
#     matchcols <- match(names(ms), names(obs))
#     if (any(!is.na(matchcols))) {
#       ms[!is.na(matchcols)] <- ms[!is.na(matchcols)] - obs[matchcols[!is.na(matchcols)]]
#       statsmatrix[,!is.na(matchcols)] <- sweep(as.matrix(
#          statsmatrix[,!is.na(matchcols)]), 2, ms[!is.na(matchcols)], "-")
#     }
#   }
  list(statsmatrix=statsmatrix, newnetwork=newnetwork, 
       meanstats=Clist$meanstats, nedges=network.edgecount(newnetwork))
}





###############################################################################
# The <ergm.mcmcslave> function is that which the slaves will call to perform
# a validation on the mcmc equal to their slave number. It also returns an
# MCMC sample.
#
# --PARAMETERS--
#   Clist     : the list of parameters returned by <ergm.Cprepare>
#   MHproposal: the MHproposal list as returned by <getMHproposal>
#   eta0      : the canonical eta parameters
#   MCMCparams: a list of parameters for controlling the MCMC algorithm;
#               recognized components include:
#       samplesize  : the number of networks to be sampled
#       stats       : ??
#       interval    : the number of proposals to ignore between sampled networks
#       burnin      : the number of proposals to initially ignore for the burn-in
#                     period
#   maxedges  : the maximum number of edges?? - this is merely used to indicate
#               whether the new generated network should be recorded in the C
#               code to pass back to R;  non-negative integers imply 'yes', 0
#               implies 'no'
#   verbose   : whether the C code should be verbose (T or F) 
#
# --RETURNED--
#   the MCMC sample as a list of the following:
#     s         : the statsmatrix
#     newnwtails: the vector of tails for the new network- is this the final
#                 network sampled? - is this the original nw if 'maxedges' is 0
#     newnwheads: the vector of heads for the new network - same q's
#
###############################################################################

ergm.mcmcslave <- function(Clist,MHproposal,eta0,MCMCparams,maxedges,verbose) {
  # A subroutine to allow caller to override some settings or resume
  # from a pervious run.
  dorun <- function(prev.run=NULL, burnin=NULL, samplesize=NULL, interval=NULL){
    
    numnetworks <- 0

    if(is.null(prev.run)){ # Start from Clist
      nedges <- c(Clist$nedges,0,0)
      tails <- Clist$tails
      heads <- Clist$heads
      weights <- Clist$weights
      stats <- rep(0, Clist$nstats)
    }else{ # Pick up where we left off
      nedges <- prev.run$newnwtails[1]
      tails <- prev.run$newnwtails[2:(nedges+1)]
      heads <- prev.run$newnwheads[2:(nedges+1)]
      weights <- prev.run$newnwweights[2:(nedges+1)]
      nedges <- c(nedges,0,0)
      stats <- matrix(prev.run$s,
                      ncol=Clist$nstats,
                      byrow = TRUE)
      stats <- stats[nrow(stats),]
    }

    if(is.null(burnin)) burnin <- MCMCparams$burnin
    if(is.null(samplesize)) samplesize <- MCMCparams$samplesize
    if(is.null(interval)) interval <- MCMCparams$interval
    
    if(is.null(Clist$weights)){
      z <- .C("MCMC_wrapper",
              as.integer(numnetworks), as.integer(nedges),
              as.integer(tails), as.integer(heads),
              as.integer(Clist$maxpossibleedges), as.integer(Clist$n),
              as.integer(Clist$dir), as.integer(Clist$bipartite),
              as.integer(Clist$nterms),
              as.character(Clist$fnamestring),
              as.character(Clist$snamestring),
              as.character(MHproposal$name), as.character(MHproposal$package),
              as.double(c(Clist$inputs,MHproposal$inputs)), as.double(eta0),
              as.integer(samplesize),
              s = as.double(rep(stats, samplesize)),
              as.integer(burnin), 
              as.integer(interval),
              newnwtails = integer(maxedges),
              newnwheads = integer(maxedges),
              as.integer(verbose), as.integer(MHproposal$arguments$constraints$bd$attribs),
              as.integer(MHproposal$arguments$constraints$bd$maxout), as.integer(MHproposal$arguments$constraints$bd$maxin),
              as.integer(MHproposal$arguments$constraints$bd$minout), as.integer(MHproposal$arguments$constraints$bd$minin),
              as.integer(MHproposal$arguments$constraints$bd$condAllDegExact), as.integer(length(MHproposal$arguments$constraints$bd$attribs)),
              as.integer(maxedges),
              status = integer(1),
              PACKAGE="ergm")
      
      # save the results
      list(s=z$s, newnwtails=z$newnwtails, newnwheads=z$newnwheads, status=z$status)
    }else{
      z <- .C("WtMCMC_wrapper",
              as.integer(length(nedges)), as.integer(nedges),
              as.integer(tails), as.integer(heads), as.double(weights),
              as.integer(Clist$maxpossibleedges), as.integer(Clist$n),
              as.integer(Clist$dir), as.integer(Clist$bipartite),
              as.integer(Clist$nterms),
              as.character(Clist$fnamestring),
              as.character(Clist$snamestring),
              as.character(MHproposal$name), as.character(MHproposal$package),
              as.double(c(Clist$inputs,MHproposal$inputs)), as.double(eta0),
              as.integer(samplesize),
              s = as.double(rep(stats, samplesize)),
              as.integer(burnin), 
              as.integer(interval),
              newnwtails = integer(maxedges),
              newnwheads = integer(maxedges),
              newnwweights = double(maxedges),
              as.integer(verbose), 
              as.integer(maxedges),
              status = integer(1),
              PACKAGE="ergm")
      # save the results
      list(s=z$s, newnwtails=z$newnwtails, newnwheads=z$newnwheads, newnwweights=z$newnwweights, status=z$status)
    }
  }
  
  if(!is.null(MCMCparams$burnin.retries) && MCMCparams$burnin.retries>0){
    library(coda)

    out <- NULL

    for(try in seq_len(MCMCparams$burnin.retries+1)){
      samplesize <- min(MCMCparams$samplesize,MCMCparams$burnin*MCMCparams$burnin.check.last)
      burnin <- ceiling(MCMCparams$burnin*(1-MCMCparams$burnin.check.last))
      interval <- ceiling(MCMCparams$burnin*MCMCparams$burnin.check.last/samplesize)
      out<-dorun(prev.run=out,
                 burnin = burnin,
                 samplesize = samplesize,
                 interval = interval)
      # Stop if something went wrong.
      if(out$status!=0) return(out)

      # Get the vector of the last burnin draws.
      burnin.stats <- matrix(out$s, nrow=samplesize,
                             ncol=Clist$nstats,
                             byrow = TRUE)[,Clist$diagnosable,drop=FALSE]
      colnames(burnin.stats) <- names(Clist$diagnosable)[Clist$diagnosable==TRUE]

      if(MCMCparams$runtime.traceplot) plot(mcmc(burnin.stats,start=burnin+1,burnin+samplesize*interval,thin=interval),ask=FALSE,smooth=TRUE,density=FALSE)
      
      # Coda's implementation uses spectrum0, which is not robust enough.
      my.geweke.diag<-function (x, frac1 = 0.1, frac2 = 0.5){
        x <- as.mcmc(x)
        xstart <- c(start(x), end(x) - frac2 * (end(x) - start(x)))
        xend <- c(start(x) + frac1 * (end(x) - start(x)), end(x))
        y.variance <- y.mean <- vector("list", 2)
        for (i in 1:2) {
          y <- window(x, start = xstart[i], end = xend[i])
          y.mean[[i]] <- apply(as.matrix(y), 2, mean)
          y.variance[[i]] <- spectrum0.ar(y)$spec/niter(y)
        }
        z <- (y.mean[[1]] - y.mean[[2]])/sqrt(y.variance[[1]] + y.variance[[2]])

        # Output 2-sided P-value, rather than Z score.
        pnorm(abs(z),0,1,lower.tail=FALSE)*2
               
      }

      # Bonferroni adjustment
      failed <- my.geweke.diag(burnin.stats) < MCMCparams$burnin.check.alpha/ncol(burnin.stats)
      # If a statistic didn't mix at all, fail it as well.
      failed[is.na(failed)] <- TRUE
      if(any(failed)){
        cat("Burn-in failed to converge or mixed very poorly for statistics", paste(names(Clist$diagnosable[Clist$diagnosable])[failed],collapse=", "), ". Rerunning.\n")
        if(try == MCMCparams$burnin.retries+1) burnin.failed <- TRUE
      }
      else{
        burnin.failed <- FALSE
        break
      }
    }

    # Do the actual sampling run. Note that we've already done the burnin.
    out <- dorun(prev.run=out, burnin=0)
    if(MCMCparams$runtime.traceplot) {
      stats <- matrix(out$s, nrow=MCMCparams$samplesize,
                      ncol=Clist$nstats,
                      byrow = TRUE)[,Clist$diagnosable,drop=FALSE]
      colnames(stats) <- names(Clist$diagnosable)[Clist$diagnosable==TRUE]
      
      plot(mcmc(stats,start=1,end=MCMCparams$samplesize*MCMCparams$interval,thin=MCMCparams$interval),ask=FALSE,smooth=TRUE,density=FALSE)
    }
    out$burnin.failed <- burnin.failed
    out
  } else {
    out <- dorun()
    if(MCMCparams$runtime.traceplot) {
      stats <- matrix(out$s, nrow=MCMCparams$samplesize,
                      ncol=Clist$nstats,
                      byrow = TRUE)[,Clist$diagnosable,drop=FALSE]
      colnames(stats) <- names(Clist$diagnosable)[Clist$diagnosable==TRUE]
      
      plot(mcmc(stats,start=MCMCparams$burnin+1,MCMCparams$burnin+MCMCparams$samplesize*MCMCparams$interval,thin=MCMCparams$interval),ask=FALSE,smooth=TRUE,density=FALSE)
    }
    out
  }
}
