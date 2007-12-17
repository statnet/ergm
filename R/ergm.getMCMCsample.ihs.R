ergm.getMCMCsample.ihs <- ergm.getMCMCsample <- function(nw, model, MHproposal, eta0, MCMCparams, 
                               verbose) {
# Note:  In reality, there should be many fewer arguments to this function,
# since most info should be passed via Clist (this is, after all, what Clist
# is for:  Holding all arguments required for the .C call).  In particular,
# the elements of MHproposal, MCMCparams, verbose should certainly
# be part of Clist.  But this is a project for another day!
  Clist <- ergm.Cprepare(nw, model)
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
    z <- .C("MCMC_wrapper",
            as.integer(Clist$heads), as.integer(Clist$tails), 
            as.integer(Clist$nedges), as.integer(Clist$n),
            as.integer(Clist$dir), as.integer(Clist$bipartite),
            as.integer(Clist$nterms), 
            as.character(Clist$fnamestring),
            as.character(Clist$snamestring),
            as.character(MHproposal$name), as.character(MHproposal$package),
#  Add:  as.double(length(MHproposal$args)), as.double(MHproposal$args), 
            as.double(Clist$inputs), as.double(eta0),
            as.integer(MCMCparams$samplesize),
            s = as.double(t(MCMCparams$stats)),
            as.integer(MCMCparams$burnin), as.integer(MCMCparams$interval), 
            newnwheads = integer(maxedges),
            newnwtails = integer(maxedges), 
            as.integer(verbose), as.integer(MHproposal$bd$attribs), 
            as.integer(MHproposal$bd$maxout), as.integer(MHproposal$bd$maxin),
            as.integer(MHproposal$bd$minout), as.integer(MHproposal$bd$minin),
            as.integer(MHproposal$bd$condAllDegExact), as.integer(length(MHproposal$bd$attribs)), 
            as.integer(maxedges),
            as.integer(MCMCparams$Clist.miss$heads), as.integer(MCMCparams$Clist.miss$tails),
            as.integer(MCMCparams$Clist.miss$nedges),
            PACKAGE="ergm") 
    if(z$newnwheads[1] >= 50000-1){
      stop(paste("\n The network has more than 50000 edges, and the model is likely to be degenerate.\n",
                  "Try starting the algorithm at an alternative model\n",
                  "(That is, changing the model terms or the 'theta0' argument).\n",
                  "The current theta0 is:\n"))
                  print(theta0)
    }
    statsmatrix <- matrix(z$s, nrow=MCMCparams$samplesize,
                          ncol=Clist$nparam,
                          byrow = TRUE)

    newnetwork <- newnw.extract(nw,z)

  }else{
    if(verbose){cat("Engaging warp drive ...\n")}
    MCMCparams.parallel <- MCMCparams
    MCMCparams.parallel$samplesize <- round(MCMCparams$samplesize / MCMCparams$parallel)
    MCMCparams.parallel$stats <- MCMCparams$stats[1:MCMCparams.parallel$samplesize,]
    require(snow)
#
# Start PVM if necessary
#
    if(getClusterOption("type")=="PVM"){
     require(rpvm)
     PVM.running <- try(.PVM.config(), silent=TRUE)
     if(inherits(PVM.running,"try-error")){
      hostfile <- paste(Sys.getenv("HOME"),"/.xpvm_hosts",sep="")
      .PVM.start.pvmd(hostfile)
      cat("no problem... PVM started by R...\n")
     }
    }
#
#   Start Cluster
#
    cl<-makeCluster(MCMCparams$parallel)
    clusterSetupRNG(cl)
    clusterEvalQ(cl,library(ergm))
#   clusterEvalQ(cl,eval(paste("library(",packagename,")",sep="")))
#
# Start PVM if necessary
#
    if(getClusterOption("type")=="PVM"){
     PVM.running <- try(.PVM.config(), silent=TRUE)
     if(inherits(PVM.running,"try-error")){
      hostfile <- paste(Sys.getenv("HOME"),"/.xpvm_hosts",sep="")
      .PVM.start.pvmd(hostfile)
      cat("no problem... PVM started by R...\n")
     }
    }
#
#   Run the jobs with rpvm or Rmpi
#
    outlist <- clusterCall(cl,mcmcslave,
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
       ncol=Clist$nparam,
       byrow = TRUE))
#    if(z$newnw[1]>1){
#      newedgelist <- rbind(newedgelist,
     #                           matrix(z$newnw[2:z$newnw[1]], ncol=2, byrow=TRUE))
   }
    newnetwork<-newnw.extract(nw,z)
    cat("parallel samplesize=",nrow(statsmatrix),"by",
        MCMCparams.parallel$samplesize,"\n")
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
  list(statsmatrix=statsmatrix, newnetwork=newnetwork, meanstats=Clist$meanstats)
}
