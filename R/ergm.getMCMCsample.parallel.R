ergm.getMCMCsample.parallel <- function(nw, model, MHproposal, eta0, MCMCparams, 
                               verbose, response=NULL) {
# Note:  In reality, there should be many fewer arguments to this function,
# since most info should be passed via Clist (this is, after all, what Clist
# is for:  Holding all arguments required for the .C call).  In particular,
# the elements of MHproposal, MCMCparams, verbose should certainly
# be part of Clist.  But this is a project for another day!
  Clist <- ergm.Cprepare(nw, model, response=response)
  if(is.null(MCMCparams$Clist.dt)){
    MCMCparams$Clist.dt <- list(heads=NULL, tails=NULL, nedges=0, dir=is.directed(nw))
  }
  maxedges <- max(5000, Clist$nedges)
#
#   Check for truncation of the returned edge list
#
  z <- list(newnwheads=maxedges+1)
  while(z$newnwheads[1] >= maxedges){
   maxedges <- 10*maxedges
#
#  Parallel running
#
   if(MCMCparams$parallel==0){
    flush.console()
    z <- ergm.mcmcslave(Clist,MHproposal,eta0,MCMCparams,maxedges,verbose)
    nedges <- z$newnwheads[1]
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
    nedges <- z$newnwheads[1]
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
# Function the slaves will call to perform a validation on the
# mcmc equal to their slave number.
# Assumes: Clist MHproposal eta0 MCMCparams maxedges verbose
ergm.mcmcslave <- function(Clist,MHproposal,eta0,MCMCparams,maxedges,verbose) {
  numnetworks <- 0
  nedges <- c(Clist$nedges,0,0)
  heads <- Clist$heads
  tails <- Clist$tails
  if(!is.null(MCMCparams$Clist.miss)){
    nedges[2] <- MCMCparams$Clist.miss$nedges
    heads <- c(heads, MCMCparams$Clist.miss$heads)
    tails <- c(tails, MCMCparams$Clist.miss$tails)
  }
  if(!is.null(MCMCparams$Clist.dt)){
    nedges[3] <- MCMCparams$Clist.dt$nedges
    heads <- c(heads, MCMCparams$Clist.dt$heads)
    tails <- c(tails, MCMCparams$Clist.dt$tails)
  }
  if(is.null(Clist$weights)){
  z <- .C("MCMC_wrapper",
  as.integer(numnetworks), as.integer(nedges),
  as.integer(heads), as.integer(tails),
  as.integer(Clist$maxpossibleedges), as.integer(Clist$n),
  as.integer(Clist$dir), as.integer(Clist$bipartite),
  as.integer(Clist$nterms),
  as.character(Clist$fnamestring),
  as.character(Clist$snamestring),
  as.character(MHproposal$name), as.character(MHproposal$package),
  as.double(Clist$inputs), as.double(eta0),
  as.integer(MCMCparams$samplesize),
  s = as.double(t(MCMCparams$stats)),
  as.integer(MCMCparams$burnin), 
  as.integer(MCMCparams$interval),
  newnwheads = integer(maxedges),
  newnwtails = integer(maxedges),
  as.integer(verbose), as.integer(MHproposal$bd$attribs),
  as.integer(MHproposal$bd$maxout), as.integer(MHproposal$bd$maxin),
  as.integer(MHproposal$bd$minout), as.integer(MHproposal$bd$minin),
  as.integer(MHproposal$bd$condAllDegExact), as.integer(length(MHproposal$bd$attribs)),
  as.integer(maxedges),
  PACKAGE="ergm")

  # save the results
    list(s=z$s, newnwheads=z$newnwheads, newnwtails=z$newnwtails)
  }else{
    z <- .C("WtMCMC_wrapper",
            as.integer(numnetworks), as.integer(nedges),
            as.integer(Clist$heads), as.integer(Clist$tails), as.double(Clist$weights),
            as.integer(Clist$maxpossibleedges), as.integer(Clist$n),
            as.integer(Clist$dir), as.integer(Clist$bipartite),
            as.integer(Clist$nterms),
            as.character(Clist$fnamestring),
            as.character(Clist$snamestring),
            as.character(MHproposal$name), as.character(MHproposal$package),
            as.double(Clist$inputs), as.double(eta0),
            as.integer(MCMCparams$samplesize),
            s = as.double(t(MCMCparams$stats)),
            as.integer(MCMCparams$burnin), 
            as.integer(MCMCparams$interval),
            newnwheads = integer(maxedges),
            newnwtails = integer(maxedges),
            newnwweights = double(maxedges),
            as.integer(verbose), 
            as.integer(maxedges),
            PACKAGE="ergm")
    # save the results
    list(s=z$s, newnwheads=z$newnwheads, newnwtails=z$newnwtails, newnwweights=z$newnwweights)
  }
}
