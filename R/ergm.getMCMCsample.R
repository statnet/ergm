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
                                        verbose, response=NULL, ...) {
  
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
    z <- ergm.mcmcslave(Clist,MHproposal,eta0,control,verbose,...)
    
    if(z$status == 1){ # MCMC_TOO_MANY_EDGES, exceeding even control$MCMC.max.maxedges
      return(list(status=z$status))
    }
    
    if(z$status == 2){ # MCMC_MH_FAILED
      # MH proposal failed somewhere. Throw an error.
      stop("Sampling failed due to a Metropolis-Hastings proposal failing.")
    }
    
    if(!is.null(z$burnin.failed) && z$burnin.failed) warning("Burn-in failed to converge after retries.")
    
    statsmatrix <- matrix(z$s,
                          ncol=Clist$nstats,
                          byrow = TRUE)
    newnetwork <- newnw.extract(nw,z,response=response,
                                output=control$network.output)
    newnetworks <- list(newnetwork)

    final.interval <- z$final.interval
  }else{
    control.parallel <- control
    control.parallel$MCMC.samplesize <- ceiling(control$MCMC.samplesize / control$parallel)
    control.parallel$MCMC.effectiveSize <- ceiling(control$MCMC.effectiveSize / control$parallel)
    
    cl <- ergm.getCluster(control, verbose)
    #
    #   Run the jobs with rpvm or Rmpi
    #
    flush.console()
    outlist <- {
      if(nnw>1) clusterMap(cl,ergm.mcmcslave,
                           Clist, MoreArgs=list(MHproposal,eta0,control.parallel,verbose,...))
      else clusterCall(cl,ergm.mcmcslave,
                       Clist,MHproposal,eta0,control.parallel,verbose,...)
    }
    #
    #   Process the results
    #
    statsmatrix <- NULL
    newnetworks <- list()
    final.interval <- c()
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
                           matrix(z$s,
                                  ncol=(if(nnw==1) Clist else Clist[[i]])$nstats,
                                  byrow = TRUE))
      newnetworks[[i]]<-newnw.extract(if(nnw==1) nw else nw[[i]],z,
                                  response=response,output=control$network.output)
      final.interval <- c(final.interval, z$final.interval)
    }
    newnetwork<-newnetworks[[1]]

    if(verbose){cat("parallel samplesize=",nrow(statsmatrix),"by",
                    control.parallel$MCMC.samplesize,"\n")}
    
    ergm.stopCluster(cl)
  }

  colnames(statsmatrix) <- model$coef.names

  statsmatrix[is.na(statsmatrix)] <- 0
  list(statsmatrix=statsmatrix, newnetwork=newnetwork, newnetworks=newnetworks, status=0, final.interval=final.interval)
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
#   ...       : optional arguments
#
# --RETURNED--
#   the MCMC sample as a list of the following:
#     s         : the statsmatrix
#     newnwtails: the vector of tails for the new network- is this the final
#                 network sampled? - is this the original nw if 'maxedges' is 0
#     newnwheads: the vector of heads for the new network - same q's
#
###############################################################################

ergm.mcmcslave <- function(Clist,MHproposal,eta0,control,verbose,...) {
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
                as.double(c(Clist$inputs,MHproposal$inputs)), as.double(.deinf(eta0)),
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
        
        # save the results (note that if prev.run is NULL, c(NULL$s,z$s)==z$s.
        z<-list(s=c(prev.run$s,z$s), newnwtails=z$newnwtails, newnwheads=z$newnwheads, status=z$status, maxedges=maxedges)
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
                as.double(c(Clist$inputs,MHproposal$inputs)), as.double(.deinf(eta0)),
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
        z<-list(s=c(prev.run$s,z$s), newnwtails=z$newnwtails, newnwheads=z$newnwheads, newnwweights=z$newnwweights, status=z$status, maxedges=maxedges)
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

  if(!is.null(control$MCMC.effectiveSize)){
    if(verbose) cat("Beginning adaptive MCMC...\n")

    
    interval <- control$MCMC.interval
    stats <- out <- NULL
    for(mcrun in seq_len(control$MCMC.effectiveSize.maxruns)){
      out<-dorun(prev.run=out,
                 burnin = interval, # I.e., skip that much before the first draw.
                 samplesize = control$MCMC.effectiveSize,
                 interval = interval,
                 maxedges = out$maxedges
                 )
      
      # Stop if something went wrong.
      if(out$status!=0) return(out)
      
      stats <- rbind(matrix(out$s, # s is cumulative
                            ncol=Clist$nstats,
                            byrow = TRUE)
                     )
      colnames(stats) <- names(Clist$diagnosable)
      
      if(nrow(stats)>=control$MCMC.effectiveSize*8){
        stats <- stats[seq_len(floor(nrow(stats)/2))*2,,drop=FALSE]
        out$s <- c(t(stats))
        interval <- interval*2
        if(verbose) cat("Increasing thinning to",interval,".\n")
      }
      
      esteq <-
        if(all(c("theta","etamap") %in% names(list(...)))) .ergm.esteq(list(...)$theta, list(etamap=list(...)$etamap), stats)
        else stats[,Clist$diagnosable,drop=FALSE]
      
      meS <- .max.effectiveSize(esteq)
      if(verbose) cat("Maximum Harmonic mean ESS of",meS$eS,"attained with burn-in of", round(meS$b/nrow(stats)*100),"%.\n")
      
      if(meS$eS>=control$MCMC.effectiveSize){
        if(verbose) cat("Target ESS achieved. Returning.\n")

        if(meS$burnin) stats <- stats[-seq_len(meS$burnin),,drop=FALSE]
        out$s <- c(t(stats))
        out$final.interval <- interval
        return(out)
      }
      
      if(control$MCMC.runtime.traceplot) plot(window(mcmc(esteq,thin=max(1,floor(nrow(esteq)/1000)))),ask=FALSE,smooth=TRUE,density=FALSE)
      
    }
    warning("Unable to reach target effective size in iterations alotted.")
    
    out$final.interval <- interval
    if(meS$burnin) stats <- stats[-seq_len(meS$burnin),,drop=FALSE]
    out$s <- c(t(stats))
    out
  }else dorun()  
}


.max.effectiveSize <- function(x, npts=20, base=3/4){
  es <- function(b){
    if(b>0) x <- x[-seq_len(b),,drop=FALSE]
    effSizes <- effectiveSize(x)
    mean.fn <- function(x) x^(-1)
    mean.ifn <- function(x) x^(-1)
    mean.ifn(mean(mean.fn(effSizes)))
  }

  pts <- sort(round(base^seq_len(npts)*nrow(x)))
  ess <- sapply(pts, es)

  list(burnin=pts[which.max(ess)], eS=max(ess))
}
