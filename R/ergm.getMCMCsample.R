#  File R/ergm.getMCMCsample.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
#==============================================================================
# This file contains the 2 following functions for getting an MCMC sample
#      <ergm.getMCMCsample>
#      <ergm.mcmcslave>
#==============================================================================




#########################################################################################
# The <ergm.getMCMCsample> function samples networks using an MCMC algorithm via
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
#   control:  list of MCMC tuning parameters; those recognized include
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
# the elements of MHproposal, control, verbose should certainly
# be part of Clist.  But this is a project for another day!
#
# --RETURNED--
#   the sample as a list containing:
#     statsmatrix:  the stats matrix for the sampled networks, RELATIVE TO THE ORIGINAL
#                   NETWORK!
#     newnetwork :  the edgelist of the final sampled network
#     nedges     :  the number of edges in the 'newnetwork' ??
#
#########################################################################################

ergm.getMCMCsample <- function(nw, model, MHproposal, eta0, control, 
                                        verbose, response=NULL) {
  
  if(is.network(nw[[1]])){ # I.e., we are dealing with a list of initial networks.
    nnw <- length(nw)
    # FIXME: This doesn't have to be the case:
    if(control$parallel!=nnw)
      stop("Number of initial networks passed to ergm.getMCMCsample must equal control$parallel (for now).")
  }else nnw <- 1
     
  Clist <- if(nnw>1) lapply(nw, ergm.Cprepare, model,
response=response) else ergm.Cprepare(nw, model, response=response)

  if(control$parallel==0){
    flush.console()
    z <- ergm.mcmcslave(Clist,MHproposal,eta0,control,verbose)
    
    if(z$status == 1){ # MCMC_TOO_MANY_EDGES, exceeding even control$MCMC.max.maxedges
      return(list(status=z$status))
    }
    
    if(z$status == 2){ # MCMC_MH_FAILED
      # MH proposal failed somewhere. Throw an error.
      stop("Sampling failed due to a Metropolis-Hastings proposal failing.")
    }
    
    if(!is.null(z$burnin.failed) && z$burnin.failed) warning("Burn-in failed to converge after retries.")
    
    statsmatrix <- matrix(z$s, nrow=control$MCMC.samplesize,
                          ncol=Clist$nstats,
                          byrow = TRUE)
    newnetwork <- newnw.extract(nw,z,response=response,
                                output=control$network.output)
    newnetworks <- list(newnetwork)
  }else{
    control.parallel <- control
    control.parallel$MCMC.samplesize <- ceiling(control$MCMC.samplesize / control$parallel)
    
    cl <- ergm.getCluster(control, verbose)
    #
    #   Run the jobs with rpvm or Rmpi
    #
    flush.console()
    outlist <- {
      if(nnw>1) clusterMap(cl,ergm.mcmcslave,
                           Clist, MoreArgs=list(MHproposal,eta0,control.parallel,verbose))
      else clusterCall(cl,ergm.mcmcslave,
                       Clist,MHproposal,eta0,control.parallel,verbose)
    }
    #
    #   Process the results
    #
    statsmatrix <- NULL
    newnetworks <- list()
    for(i in (1:control$parallel)){
      z <- outlist[[i]]

      if(z$status == 1){ # MCMC_TOO_MANY_EDGES, exceeding even control$MCMC.max.maxedges
        return(list(status=z$status))
      }
      
      if(z$status == 2){ # MCMC_MH_FAILED
        # MH proposal failed somewhere. Throw an error.
        stop("Sampling failed due to a Metropolis-Hastings proposal failing.")
      }
      
      if(!is.null(z$burnin.failed) && z$burnin.failed) warning("Burn-in failed to converge after retries.")
      
      statsmatrix <- rbind(statsmatrix,
                           matrix(z$s, nrow=control.parallel$MCMC.samplesize,
                                  ncol=(if(nnw==1) Clist else Clist[[i]])$nstats,
                                  byrow = TRUE))
      newnetworks[[i]]<-newnw.extract(if(nnw==1) nw else nw[[i]],z,response=response)
    }
    newnetwork<-newnetworks[[1]]

    if(verbose){cat("parallel samplesize=",nrow(statsmatrix),"by",
                    control.parallel$MCMC.samplesize,"\n")}
    
    ergm.stopCluster(cl)
  }

  colnames(statsmatrix) <- model$coef.names

  statsmatrix[is.na(statsmatrix)] <- 0
  list(statsmatrix=statsmatrix, newnetwork=newnetwork, newnetworks=newnetworks, status=0)
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
#   control: a list of parameters for controlling the MCMC algorithm;
#               recognized components include:
#       samplesize  : the number of networks to be sampled
#       interval    : the number of proposals to ignore between sampled networks
#       burnin      : the number of proposals to initially ignore for the burn-in
#                     period
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

ergm.mcmcslave <- function(Clist,MHproposal,eta0,control,verbose) {
  # A subroutine to allow caller to override some settings or resume
  # from a pervious run.
  dorun <- function(prev.run=NULL, burnin=NULL, samplesize=NULL, interval=NULL, maxedges=NULL){
    
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

    if(is.null(burnin)) burnin <- control$MCMC.burnin
    if(is.null(samplesize)) samplesize <- control$MCMC.samplesize
    if(is.null(interval)) interval <- control$MCMC.interval
    if(is.null(maxedges)) maxedges <- control$MCMC.init.maxedges

    repeat{
      if(is.null(Clist$weights)){
        z <- .C("MCMC_wrapper",
                as.integer(numnetworks), as.integer(nedges),
                as.integer(tails), as.integer(heads),
                as.integer(Clist$n),
                as.integer(Clist$dir), as.integer(Clist$bipartite),
                as.integer(Clist$nterms),
                as.character(Clist$fnamestring),
                as.character(Clist$snamestring),
                as.character(MHproposal$name), as.character(MHproposal$pkgname),
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
        z<-list(s=z$s, newnwtails=z$newnwtails, newnwheads=z$newnwheads, status=z$status, maxedges=maxedges)
      }else{
        z <- .C("WtMCMC_wrapper",
                as.integer(length(nedges)), as.integer(nedges),
                as.integer(tails), as.integer(heads), as.double(weights),
                as.integer(Clist$n),
                as.integer(Clist$dir), as.integer(Clist$bipartite),
                as.integer(Clist$nterms),
                as.character(Clist$fnamestring),
                as.character(Clist$snamestring),
                as.character(MHproposal$name), as.character(MHproposal$pkgname),
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
        z<-list(s=z$s, newnwtails=z$newnwtails, newnwheads=z$newnwheads, newnwweights=z$newnwweights, status=z$status, maxedges=maxedges)
      }
      if(z$status!=1) return(z) # Handle everything except for MCMC_TOO_MANY_EDGES elsewhere.

      # The following is only executed (and the loop continued) if too many edges.
      maxedges <- maxedges * 10
      if(!is.null(control$MCMC.max.maxedges)){
        if(maxedges == control$MCMC.max.maxedges*10) # True iff the previous maxedges exactly equaled control$MCMC.max.maxedges and that was too small.
          return(z) # This will kick the too many edges problem upstairs, so to speak.
        maxedges <- min(maxedges, control$MCMC.max.maxedges)
      }

    }
  }
  
  if(!is.null(control$MCMC.burnin.retries) && control$MCMC.burnin.retries>0){
    out <- NULL

    for(try in seq_len(control$MCMC.burnin.retries+1)){
      samplesize <- min(control$MCMC.samplesize,control$MCMC.burnin*control$MCMC.burnin.check.last)
      burnin <- ceiling(control$MCMC.burnin*(1-control$MCMC.burnin.check.last))
      interval <- ceiling(control$MCMC.burnin*control$MCMC.burnin.check.last/samplesize)
      out<-dorun(prev.run=out,
                 burnin = burnin,
                 samplesize = samplesize,
                 interval = interval,
                 maxedges = out$maxedges # note that the first time through, maxedges=NULL$maxedges, which is NULL.
                 )
      # Stop if something went wrong.
      if(out$status!=0) return(out)

      # Get the vector of the last burnin draws.
      burnin.stats <- matrix(out$s, nrow=samplesize,
                             ncol=Clist$nstats,
                             byrow = TRUE)[,Clist$diagnosable,drop=FALSE]
      colnames(burnin.stats) <- names(Clist$diagnosable)[Clist$diagnosable==TRUE]

      if(control$MCMC.runtime.traceplot) plot(mcmc(burnin.stats,start=burnin+1,burnin+samplesize*interval,thin=interval),ask=FALSE,smooth=TRUE,density=FALSE)
      
      # Bonferroni adjustment
      failed <- geweke.diag.ar(burnin.stats)$p.val < control$MCMC.burnin.check.alpha/ncol(burnin.stats)
      # If a statistic didn't mix at all, fail it as well.
      failed[is.na(failed)] <- TRUE
      if(any(failed)){
        cat("Burn-in failed to converge or mixed very poorly for statistics", paste.and(names(Clist$diagnosable[Clist$diagnosable])[failed]), ". Rerunning.\n")
        if(try == control$MCMC.burnin.retries+1) burnin.failed <- TRUE
      }
      else{
        burnin.failed <- FALSE
        break
      }
    }

    # Do the actual sampling run. Note that we've already done the burnin.
    out <- dorun(prev.run=out, burnin=0, maxedges=out$maxedges)
    if(control$MCMC.runtime.traceplot) {
      stats <- matrix(out$s, nrow=control$MCMC.samplesize,
                      ncol=Clist$nstats,
                      byrow = TRUE)[,Clist$diagnosable,drop=FALSE]
      colnames(stats) <- names(Clist$diagnosable)[Clist$diagnosable==TRUE]
      
      plot(mcmc(stats,start=1,end=control$MCMC.samplesize*control$MCMC.interval,thin=control$MCMC.interval),ask=FALSE,smooth=TRUE,density=FALSE)
    }
    out$burnin.failed <- burnin.failed
    out
  } else {
    out <- dorun()
    if(control$MCMC.runtime.traceplot) {
      stats <- matrix(out$s, nrow=control$MCMC.samplesize,
                      ncol=Clist$nstats,
                      byrow = TRUE)[,Clist$diagnosable,drop=FALSE]
      colnames(stats) <- names(Clist$diagnosable)[Clist$diagnosable==TRUE]
      
      plot(mcmc(stats,start=control$MCMC.burnin+1,control$MCMC.burnin+control$MCMC.samplesize*control$MCMC.interval,thin=control$MCMC.interval),ask=FALSE,smooth=TRUE,density=FALSE)
    }
    out
  }
}
