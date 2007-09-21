ergm.phase12.dyn <- function(g, model.form, model.diss, 
                        MHproposal.form, MHproposal.diss, eta0, gamma0,
                        MCMCparams, verbose, BD) {
# ms <- MCMCparams$meanstats
# if(!is.null(ms)) {
#   if (is.null(names(ms)) && length(ms) == length(model.form$coef.names))
#     names(ms) <- model.form$coef.names
#   obs <- MCMCparams$orig.obs
#   obs <- obs[match(names(ms), names(obs))]
#   ms  <-  ms[match(names(obs), names(ms))]
#   matchcols <- match(names(ms), names(obs))
#   if (any(!is.na(matchcols))) {
#     ms[!is.na(matchcols)] <- ms[!is.na(matchcols)] - obs[matchcols[!is.na(matchcols)]]
#   }
# }
  Clist.form <- ergm.Cprepare(g, model.form)
  Clist.diss <- ergm.Cprepare(g, model.diss)
  maxchanges <- max(MCMCparams$maxchanges, Clist.form$nedges)/5
  MCMCparams$maxchanges <- MCMCparams$maxchanges/5
  z <- list(newnwhead=maxchanges+1)
  while(z$newnwhead[1]  >= maxchanges || 
        z$diffnwhead[1] >= maxchanges){
    maxchanges <- 5*maxchanges
    MCMCparams$maxchanges <- 5*MCMCparams$maxchanges
    if(verbose){cat(paste("MCMCDyn workspace is",maxchanges,"\n"))}
#
#  Parallel running
#
    if(MCMCparams$parallel==0){
      z <- .C("MCMCDynPhase12",
              # Observed/starting network.
              as.integer(Clist.form$heads), as.integer(Clist.form$tails), 
              as.integer(Clist.form$nedges), as.integer(Clist.form$n),
              as.integer(Clist.form$dir), as.integer(Clist.form$bipartite),
              # Order code.
              as.integer(Clist.diss$order.code),
              # Formation terms and proposals.
              as.integer(Clist.form$nterms), as.character(Clist.form$fnamestring), as.character(Clist.form$snamestring),
              as.character(MHproposal.form$name), as.character(MHproposal.form$package),
#  Add:  as.double(length(MHproposal$args)), as.double(MHproposal$args), 
              as.double(Clist.form$inputs), eta=as.double(eta0),
              # Formation parameter fitting.
              as.double(MCMCparams$gain), as.double(MCMCparams$stats[1,]),
              as.integer(MCMCparams$phase1), as.integer(MCMCparams$nsub),
              # Dissolution terms and proposals.
              as.integer(Clist.diss$nterms), as.character(Clist.diss$fnamestring), as.character(Clist.diss$snamestring),
              as.character(MHproposal.diss$name), as.character(MHproposal.diss$package),
              as.double(Clist.diss$inputs), as.double(MCMCparams$gamma0),
              # Degree bounds.
              as.integer(BD$attribs), 
              as.integer(BD$maxout), as.integer(BD$maxin),
              as.integer(BD$minout), as.integer(BD$minin),
              as.integer(BD$condAllDegExact), as.integer(length(BD$attribs)), 
              # MCMC settings.              
              as.integer(MCMCparams$samplesize), as.integer(MCMCparams$dyninterval),
              as.integer(MCMCparams$burnin), as.integer(MCMCparams$interval),
              # Space for output.
              s.form = double(MCMCparams$samplesize * Clist.form$nparam), d.form = double(MCMCparams$samplesize * Clist.diss$nparam),
              newnwhead = integer(maxchanges), newnwtail = integer(maxchanges),
              as.double(maxchanges),
              diffnwtime = integer(maxchanges), diffnwhead = integer(maxchanges), diffnwtail = integer(maxchanges),
              # Verbosity.
              as.integer(verbose), 
          PACKAGE="ergm") 
    statsmatrix <- matrix(z$s.form, nrow=MCMCparams$samplesize,
                          ncol=Clist.form$nparam,
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
     Clist.form,
     Clist.diss,
     MHproposal.form,
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
       matrix(z$s.form, nrow=MCMCparams.parallel$samplesize,
       ncol=Clist.form$nparam,
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
  newnetwork<-newnw.extract(g,z)
#   Next create the network of differences from the origianl one
    if(z$diffnwhead[1]>1){
     diffedgelist <- cbind(z$diffnwtime[2:z$diffnwtime[1]],z$diffnwhead[2:z$diffnwhead[1]],z$diffnwtail[2:z$diffnwhead[1]])
    }else{
     diffedgelist <- matrix(0, ncol=3, nrow=0)
    }
    colnames(statsmatrix) <- model.form$coef.names
  list(statsmatrix=statsmatrix, newnetwork=newnetwork,
       meanstats=MCMCparams$meanstats,
       changed=diffedgelist,
       maxchanges=MCMCparams$maxchanges,
       eta=eta)
}
