ergm.getMCMCsample <- function(nw, model, MHproposal, eta0, MCMCparams, 
                               verbose, BD) {
# Note:  In reality, there should be many fewer arguments to this function,
# since most info should be passed via Clist (this is, after all, what Clist
# is for:  Holding all arguments required for the .C call).  In particular,
# the elements of MHproposal, MCMCparams, verbose, and BD should certainly
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
            as.integer(verbose), as.integer(BD$attribs), 
            as.integer(BD$maxout), as.integer(BD$maxin),
            as.integer(BD$minout), as.integer(BD$minin),
            as.integer(BD$condAllDegExact), as.integer(length(BD$attribs)), 
            as.integer(maxedges),
            as.integer(MCMCparams$Clist.miss$heads), as.integer(MCMCparams$Clist.miss$tails),
            as.integer(MCMCparams$Clist.miss$nedges),
            PACKAGE="statnet") 
    if(z$newnwheads[1] > 50000){
      stop(paste("The network has more then 50000 edges, and the model is likely to be degenerate.\n",
                  "Try starting the algorithm at an alternative model\n",
                  "(That is, changing the model terms or the 'theta0' argument).\n"))
    }
    statsmatrix <- matrix(z$s, nrow=MCMCparams$samplesize,
                          ncol=Clist$nparam,
                          byrow = TRUE)

    newnetwork <- newnw.extract(nw,z)

  }else{
    rpvmbasename <- paste("ergm.parallel.",Sys.getpid(),sep="")
    MCMCparams.parallel <- MCMCparams
    MCMCparams.parallel$samplesize <- round(MCMCparams$samplesize / MCMCparams$parallel)
    MCMCparams.parallel$stats <- MCMCparams$stats[1:MCMCparams.parallel$samplesize,]
    require(rpvm)
    require(MASS) # needed by rpvm
#
#   Write the slave file
#
    outsetuppvm <- ergm.rpvm.setup(rpvmbasename, verbose=verbose,
                                   packagename=packagename)
#
#   Saving the common variables
#
    save(
     Clist,
     MHproposal,
     eta0,
     MCMCparams.parallel,
     maxedges, 
     verbose,
     BD, 
     file=paste(outsetuppvm$SLAVEDIR,"/",rpvmbasename,".common.RData",sep="")
    )
#
#   Run the jobs with PVM
#
    outlist <- ergm.rpvm.run(MCMCparams$parallel, rpvmbasename)
#
#   Process the results
#
    statsmatrix <- NULL
#   newedgelist <- matrix(0, ncol=2, nrow=0)
    for(i in (1:MCMCparams$parallel)){
     load(file=
      paste(outsetuppvm$SLAVEDIR,"/",rpvmbasename,".out.",i,".RData",sep=""))
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
    ergm.rpvm.clean(rpvmbasename=rpvmbasename)
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
