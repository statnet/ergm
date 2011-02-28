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
#       Clist.dt    : this is a Clist, similar to that returned by
#                     <ergm.Cprepare>, but this is for fitting dynamic models
#       Clist.miss  : a corresponding 'Clist' for the network of missing edges,
#                     as returned by <ergm.design>
#       samplesize  : the number of networks to be sampled
#       interval    : the number of proposals to ignore between sampled networks
#       burnin      : the number of proposals to initially ignore for the burn-in
#                     period
#       stats       : ??
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
  if(is.null(MCMCparams$Clist.dt)){
    MCMCparams$Clist.dt <- list(tails=NULL, heads=NULL, nedges=0, dir=is.directed(nw))
  }
  maxedges <- max(5000, Clist$nedges)
#
#   Check for truncation of the returned edge list
#
  z <- list(newnwtails=maxedges+1)
  while(z$newnwtails[1] >= maxedges){
   maxedges <- 10*maxedges
#
#  Parallel running
#
   if(MCMCparams$parallel==0){
    flush.console()
    z <- ergm.mcmcslave(Clist,MHproposal,eta0,MCMCparams,maxedges,verbose)
    nedges <- z$newnwtails[1]
    statsmatrix <- matrix(z$s, nrow=MCMCparams$samplesize,
                          ncol=Clist$nstats,
                          byrow = TRUE)
    newnetwork <- newnw.extract(nw,z,response=response)
    if(nedges >= 50000-1){
      cat("\n Warning:")
      cat("\n   The network has more than 50000 edges, and the model is likely to be degenerate.\n")
#  THE FOLLOWING THREE LINES ARE COMMENTED OUT IN branches/2.2:
      statsmatrix <- matrix(0, nrow=MCMCparams$samplesize,
                            ncol=Clist$nstats)
      newnetwork <- network.copy(nw)
    }      
  }else{
    MCMCparams.parallel <- MCMCparams
    MCMCparams.parallel$samplesize <- round(MCMCparams$samplesize / MCMCparams$parallel)
    MCMCparams.parallel$stats <- MCMCparams$stats[1:MCMCparams.parallel$samplesize,]
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
#   if("networksis" %in% MCMCparams$packagenames){
#    clusterEvalQ(cl,library(networksis))
#   }
#    clusterEvalQ(cl,eval(paste("library(",packagename,")",sep="")))
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
     statsmatrix <- rbind(statsmatrix,
       matrix(z$s, nrow=MCMCparams.parallel$samplesize,
       ncol=Clist$nstats,
       byrow = TRUE))
#    if(z$newnw[1]>1){
#      newedgelist <- rbind(newedgelist,
     #                           matrix(z$newnw[2:z$newnw[1]], ncol=2, byrow=TRUE))
   }
    nedges <- z$newnwtails[1]
    newnetwork<-newnw.extract(nw,z,response=response)
    if(verbose){cat("parallel samplesize=",nrow(statsmatrix),"by",
	MCMCparams.parallel$samplesize,"\n")}

    stopCluster(cl)
  }
  }
  colnames(statsmatrix) <- model$coef.names

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
       meanstats=Clist$meanstats, nedges=nedges)
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
#       Clist.dt    : this is a Clist, similar to that returned by
#                     <ergm.Cprepare>, but this is for fitting dynamic models
#       Clist.miss  : a corresponding 'Clist' for the network of missing edges,
#                     as returned by <ergm.design>
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
  numnetworks <- 0
  nedges <- c(Clist$nedges,0,0)
  tails <- Clist$tails
  heads <- Clist$heads
  weights <- Clist$weights
  if(!is.null(MCMCparams$Clist.miss)){
    nedges[2] <- MCMCparams$Clist.miss$nedges
    tails <- c(tails, MCMCparams$Clist.miss$tails)
    heads <- c(heads, MCMCparams$Clist.miss$heads)
    weights <- c(weights, rep(1,nedges[1]))
  }
  if(!is.null(MCMCparams$Clist.dt)){
    nedges[3] <- MCMCparams$Clist.dt$nedges
    tails <- c(tails, MCMCparams$Clist.dt$tails)
    heads <- c(heads, MCMCparams$Clist.dt$heads)
  }
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
  as.integer(MCMCparams$samplesize),
  s = as.double(t(MCMCparams$stats)),
  as.integer(MCMCparams$burnin), 
  as.integer(MCMCparams$interval),
  newnwtails = integer(maxedges),
  newnwheads = integer(maxedges),
  as.integer(verbose), as.integer(MHproposal$bd$attribs),
  as.integer(MHproposal$bd$maxout), as.integer(MHproposal$bd$maxin),
  as.integer(MHproposal$bd$minout), as.integer(MHproposal$bd$minin),
  as.integer(MHproposal$bd$condAllDegExact), as.integer(length(MHproposal$bd$attribs)),
  as.integer(maxedges),
  PACKAGE="ergm")

  # save the results
  list(s=z$s, newnwtails=z$newnwtails, newnwheads=z$newnwheads)
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
            as.integer(MCMCparams$samplesize),
            s = as.double(t(MCMCparams$stats)),
            as.integer(MCMCparams$burnin), 
            as.integer(MCMCparams$interval),
            newnwtails = integer(maxedges),
            newnwheads = integer(maxedges),
            newnwweights = double(maxedges),
            as.integer(verbose), 
            as.integer(maxedges),
            PACKAGE="ergm")
    # save the results
    list(s=z$s, newnwtails=z$newnwtails, newnwheads=z$newnwheads, newnwweights=z$newnwweights)
  }
}
