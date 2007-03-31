ergm.phase2.dyn <- function(g, model, model.dissolve, 
                        MHproposal, eta0,
                        aDdiaginv,
                        MCMCparams, verbose, BD) {
  ms <- MCMCparams$meanstats
  if(!is.null(ms)) {
    if (is.null(names(ms)) && length(ms) == length(model$coef.names))
      names(ms) <- model$coef.names
    obs <- MCMCparams$orig.obs
    obs <- obs[match(names(ms), names(obs))]
    ms  <-  ms[match(names(obs), names(ms))]
    matchcols <- match(names(ms), names(obs))
    if (any(!is.na(matchcols))) {
      ms[!is.na(matchcols)] <- ms[!is.na(matchcols)] - obs[matchcols[!is.na(matchcols)]]
    }
  }
  Clist <- ergm.Cprepare(g, model)
  Clist.dissolve <- ergm.Cprepare(g, model.dissolve)
  maxchanges <- max(MCMCparams$maxchanges, Clist$nedges)/5
  MCMCparams$maxchanges <- MCMCparams$maxchanges/5
  z <- list(newnwhead=maxchanges+1)
  while(z$newnwhead[1]  >= maxchanges || 
        z$dissnwhead[1] >= maxchanges ||
        z$diffnwhead[1] >= maxchanges){
    maxchanges <- 5*maxchanges
    MCMCparams$maxchanges <- 5*MCMCparams$maxchanges
    if(verbose){cat(paste("MCMCDyn workspace is",maxchanges,"\n"))}
#
#  Parallel running
#
    if(MCMCparams$parallel==0){
      z <- .C("MCMCDynPhase2",
          as.integer(Clist.dissolve$order.code),
          as.double(Clist$heads), as.double(Clist$tails), 
          as.double(Clist$nedges), as.double(Clist$n),
          as.integer(Clist$dir), as.double(Clist$bipartite),
          as.integer(Clist$nterms), 
          as.character(Clist$fnamestring),
          as.character(Clist$snamestring),
          as.character(MHproposal$type), as.character(MHproposal$package),
          as.double(Clist$inputs),
          eta=as.double(eta0),
          as.double(MCMCparams$gain), as.double(ms),
          as.integer(MCMCparams$phase1),
          as.integer(MCMCparams$nsub),
          as.integer(Clist.dissolve$nterms),
          as.character(Clist.dissolve$fnamestring),
          as.character(Clist.dissolve$snamestring),
          as.double(Clist.dissolve$inputs),
          as.double(MCMCparams$samplesize),
          s = double(MCMCparams$samplesize * Clist$nparam),
          as.double(MCMCparams$burnin), as.double(MCMCparams$interval),
          newnwhead = integer(maxchanges), newnwtail = integer(maxchanges), 
          diffnwtime = integer(maxchanges), diffnwhead = integer(maxchanges), diffnwtail = integer(maxchanges), 
          dissnwtime = integer(maxchanges), dissnwhead = integer(maxchanges), dissnwtail = integer(maxchanges),
          as.integer(verbose), 
          as.double(MCMCparams$gamma), as.integer(MCMCparams$dyninterval), 
          as.integer(BD$attribs), 
          as.integer(BD$maxout), as.integer(BD$maxin),
          as.integer(BD$minout), as.integer(BD$minin),
          as.integer(BD$condAllDegExact), as.integer(length(BD$attribs)), 
          as.double(maxchanges),
          as.double(0.0), as.double(0.0), 
          as.double(0.0), as.integer(0),
          PACKAGE="statnet") 
    statsmatrix <- matrix(z$s, nrow=MCMCparams$samplesize,
                          ncol=Clist$nparam,
                          byrow = TRUE)
    eta <- z$eta
    names(eta) <- names(eta0)
  }else{
    rpvmbasename <- paste("ergm.parallel.",Sys.getpid(),sep="")
    MCMCparams.parallel <- MCMCparams
    MCMCparams.parallel$samplesize <- round(MCMCparams$samplesize / MCMCparams$parallel)
    require(rpvm)
    require(MASS) # needed by rpvm
#
#   Write the slave file
#
    outsetuppvm <- ergm.rpvm.setup.dyn(rpvmbasename, verbose=verbose,
                                       packagename=packagename)
#
#   Saving the common variables
#
    save(
     Clist,
     Clist.dissolve,
     MHproposal,
     eta0,
     MCMCparams.parallel,
     maxchanges, 
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
    cat("parallel samplesize=",nrow(statsmatrix),"by",
        MCMCparams.parallel$samplesize,"\n")
    ergm.rpvm.clean(rpvmbasename=rpvmbasename)
  }
  }
#   cat(paste("z$diffnwhead = ",maxchanges,z$diffnwhead[1],"\n"))
#   cat(paste("z$dissnwhead = ",maxchanges,z$dissnwhead[1],"\n"))
#   cat(paste("z$newwhead = ",maxchanges,z$newnwhead[1],"\n"))
    if(z$newnwhead[1]>1){
#    newedgelist <- matrix(z$newnw[2:z$newnw[1]], ncol=2, byrow=TRUE)
     newedgelist <- cbind(z$newnwhead[2:z$newnwhead[1]],z$newnwtail[2:z$newnwhead[1]])
    }else{
     newedgelist <- matrix(0, ncol=2, nrow=0)
    }
    if(z$dissnwhead[1]>1){
     dissedgelist <- cbind(z$dissnwtime[2:z$dissnwtime[1]],z$dissnwhead[2:z$dissnwhead[1]],z$dissnwtail[2:z$dissnwhead[1]])
    }else{
     dissedgelist <- matrix(0, ncol=3, nrow=0)
    }
#   Next create the network of differences from the origianl one
    if(z$diffnwhead[1]>1){
     diffedgelist <- cbind(z$diffnwtime[2:z$diffnwtime[1]],z$diffnwhead[2:z$diffnwhead[1]],z$diffnwtail[2:z$diffnwhead[1]])
    }else{
     diffedgelist <- matrix(0, ncol=3, nrow=0)
    }
    colnames(statsmatrix) <- model$coef.names
  list(statsmatrix=statsmatrix, newedgelist=newedgelist, meanstats=ms,
       changed=diffedgelist, dissolved=dissedgelist,
       maxchanges=MCMCparams$maxchanges,
       eta=eta)
}
