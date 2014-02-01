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

ergm.getCDsample <- function(nw, model, MHproposal, eta0, control, 
                             verbose, response=NULL) {
    
  Clist <- ergm.Cprepare(nw, model, response=response)

  z <- NULL

  if(control$parallel==0){
    flush.console()
    z <- ergm.cdslave(Clist,MHproposal,eta0,control,verbose)
    
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
  }else{
    control.parallel <- control
    control.parallel$MCMC.samplesize <- round(control$MCMC.samplesize / control$parallel)
    
    cl <- ergm.getCluster(control, verbose)
    #
    #   Run the jobs with rpvm or Rmpi
    #
    flush.console()
    outlist <- clusterCall(cl,ergm.cdslave,
                           Clist,MHproposal,eta0,control.parallel,verbose)
    #
    #   Process the results
    #
    statsmatrix <- NULL
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
                                  ncol=Clist$nstats,
                                  byrow = TRUE))
    }

    if(verbose){cat("parallel samplesize=",nrow(statsmatrix),"by",
                    control.parallel$MCMC.samplesize,"\n")}
    
    ergm.stopCluster(cl)
  }

  colnames(statsmatrix) <- model$coef.names

  statsmatrix[is.na(statsmatrix)] <- 0
  list(statsmatrix=statsmatrix, status=0)
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

ergm.cdslave <- function(Clist,MHproposal,eta0,control,verbose) {
    # A subroutine to allow caller to override some settings or resume
    # from a pervious run.
    numnetworks <- 0
    
    nedges <- c(Clist$nedges,0,0)
    tails <- Clist$tails
    heads <- Clist$heads
    weights <- Clist$weights
    stats <- rep(0, Clist$nstats)
    
    samplesize <- control$MCMC.samplesize
    
    if(is.null(Clist$weights)){
        z <- .C("CD_wrapper",
                as.integer(numnetworks), as.integer(nedges),
                as.integer(tails), as.integer(heads),
                as.integer(Clist$n),
                as.integer(Clist$dir), as.integer(Clist$bipartite),
                as.integer(Clist$nterms),
                as.character(Clist$fnamestring),
                as.character(Clist$snamestring),
                as.character(MHproposal$name), as.character(MHproposal$pkgname),
                as.double(c(Clist$inputs,MHproposal$inputs)), as.double(.deinf(eta0)),
                as.integer(samplesize), as.integer(control$CD.nsteps),
                s = as.double(rep(stats, samplesize)),
                as.integer(verbose), as.integer(MHproposal$arguments$constraints$bd$attribs),
                as.integer(MHproposal$arguments$constraints$bd$maxout), as.integer(MHproposal$arguments$constraints$bd$maxin),
                as.integer(MHproposal$arguments$constraints$bd$minout), as.integer(MHproposal$arguments$constraints$bd$minin),
                as.integer(MHproposal$arguments$constraints$bd$condAllDegExact), as.integer(length(MHproposal$arguments$constraints$bd$attribs)),
                status = integer(1),
                PACKAGE="ergm")
        
        # save the results
        z<-list(s=z$s, newnwtails=z$newnwtails, newnwheads=z$newnwheads, status=z$status)
      }else{
        z <- .C("WtCD_wrapper",
                as.integer(length(nedges)), as.integer(nedges),
                as.integer(tails), as.integer(heads), as.double(weights),
                as.integer(Clist$n),
                as.integer(Clist$dir), as.integer(Clist$bipartite),
                as.integer(Clist$nterms),
                as.character(Clist$fnamestring),
                as.character(Clist$snamestring),
                as.character(MHproposal$name), as.character(MHproposal$pkgname),
                as.double(c(Clist$inputs,MHproposal$inputs)), as.double(.deinf(eta0)),
                as.integer(samplesize), as.integer(control$CD.nsteps),
                s = as.double(rep(stats, samplesize)),
                as.integer(verbose), 
                status = integer(1),
                PACKAGE="ergm")
        # save the results
        z<-list(s=z$s, status=z$status)
      }

    z
    
}
