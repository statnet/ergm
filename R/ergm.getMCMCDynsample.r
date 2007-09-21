ergm.getMCMCDynsample <- function(nw, model.form, model.diss,
                                  MHproposal.form, MHproposal.diss, theta0, gamma0, MCMCparams, 
                                  verbose, BD){
# Note:  In reality, there should be many fewer arguments to this function,
# since most info should be passed via Clist (this is, after all, what Clist
# is for:  Holding all arguments required for the .C call).  In particular,
# the elements of MHproposal.form, MCMCparams, verbose, and BD should certainly
# be part of Clist.  But this is a project for another day!
#
#   Check for truncation of the returned edge list
#
  Clist.form <- ergm.Cprepare(nw, model.form)
  Clist.diss <- ergm.Cprepare(nw, model.diss)
  maxchanges <- max(MCMCparams$maxchanges, Clist.form$nedges)/5
  MCMCparams$maxchanges <- MCMCparams$maxchanges/5
  z <- list(newnwheads=maxchanges+1)
  while(z$newnwheads[1]  >= maxchanges - 10 || 
        z$diffnwheads[1] >= maxchanges - 10){
    maxchanges <- 5*maxchanges
    MCMCparams$maxchanges <- 5*MCMCparams$maxchanges
    if(verbose){cat(paste("MCMCDyn workspace is",maxchanges,"\n"))}
#
#  Parallel running
#
    if(MCMCparams$parallel==0){
      z <- .C("MCMCDyn_wrapper",
              # Observed network.
              as.integer(Clist.form$heads), as.integer(Clist.form$tails), 
              as.integer(Clist.form$nedges), as.integer(Clist.form$n),
              as.integer(Clist.form$dir), as.integer(Clist.form$bipartite),
              # Ordering of formation and dissolution.
              as.integer(Clist.diss$order.code),
              # Formation terms and proposals.
              as.integer(Clist.form$nterms), 
              as.character(Clist.form$fnamestring),
              as.character(Clist.form$snamestring),
              as.character(MHproposal.form$name), as.character(MHproposal.form$package),
              as.double(Clist.form$inputs), as.double(theta0),
              #?  Add:  as.double(length(MHproposal.form$args)), as.double(MHproposal.form$args),
              # Dissolution terms and proposals.
              as.integer(Clist.diss$nterms), 
              as.character(Clist.diss$fnamestring),
              as.character(Clist.diss$snamestring),
              as.character(MHproposal.diss$name), as.character(MHproposal.diss$package),
              as.double(Clist.diss$inputs), as.double(gamma0),
              # Degree bounds.
              as.integer(BD$attribs), 
              as.integer(BD$maxout), as.integer(BD$maxin),
              as.integer(BD$minout), as.integer(BD$minin),
              as.integer(BD$condAllDegExact), as.integer(length(BD$attribs)),
              # MCMC settings.
              as.double(MCMCparams$nsteps), as.integer(MCMCparams$dyninterval),
              as.double(MCMCparams$burnin), as.double(MCMCparams$interval),
              # Space for output.
              s.form = as.double(cbind(t(MCMCparams$stats.form),matrix(0,nrow=length(model.form$coef.names),ncol=MCMCparams$nsteps))),
              s.diss = as.double(cbind(t(MCMCparams$stats.diss),matrix(0,nrow=length(model.diss$coef.names),ncol=MCMCparams$nsteps))),
              newnwheads = integer(maxchanges), newnwtails = integer(maxchanges), 
              as.double(maxchanges),
              diffnwtime = integer(maxchanges),
              diffnwheads = integer(maxchanges),
              diffnwtails = integer(maxchanges),
              as.integer(verbose), 
              PACKAGE="ergm") 
      statsmatrix.form <- matrix(z$s.form, nrow=MCMCparams$nsteps+1,
                                 ncol=Clist.form$nparam,
                                 byrow = TRUE)[-1,,drop=FALSE]
      statsmatrix.diss <- matrix(z$s.diss, nrow=MCMCparams$nsteps+1,
                                 ncol=Clist.diss$nparam,
                                 byrow = TRUE)[-1,,drop=FALSE]
    }else{
      rpvmbasename <- paste("ergm.parallel.",Sys.getpid(),sep="")
      MCMCparams.parallel <- MCMCparams
      MCMCparams.parallel$nsteps <- round(MCMCparams$nsteps / MCMCparams$parallel)
      MCMCparams.parallel$stats <- MCMCparams$stats[1:MCMCparams.parallel$nsteps,]
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
      save(Clist.form,
           Clist.diss,
           MHproposal.form,
           MHproposal.diss,
           theta0,
           gamma0,
           MCMCparams.parallel,
           maxchanges, 
           verbose,
           BD,
           file=paste(outsetuppvm$SLAVEDIR,"/",rpvmbasename,".common.RData",sep=""))
#
#   Run the jobs with PVM
#
      outlist <- ergm.rpvm.run(MCMCparams$parallel, rpvmbasename)
#
#   Process the results
#
      statsmatrix.form <- NULL
      statsmatrix.diss <- NULL
#   newedgelist <- matrix(0, ncol=2, nrow=0)
      for(i in (1:MCMCparams$parallel)){
        load(file=paste(outsetuppvm$SLAVEDIR,"/",rpvmbasename,".out.",i,".RData",sep=""))
# NOTE: rbind is S-L-O-W for large structures, so it might be quicker to "allocate" a
# statsmatrix, and copy into it. -- PK
        statsmatrix.form <- rbind(statsmatrix.form,
                                  matrix(z$s, nrow=MCMCparams.parallel$nsteps,
                                         ncol=Clist.form$nparam,
                                         byrow = TRUE))
        statsmatrix.diss <- rbind(statsmatrix.diss,
                                      matrix(z$s, nrow=MCMCparams.parallel$nsteps,
                                             ncol=Clist.diss$nparam,
                                             byrow = TRUE))
#    if(z$newnw[1]>1){
#      newedgelist <- rbind(newedgelist,
#                           matrix(z$newnw[2:z$newnw[1]], ncol=2, byrow=TRUE))
      }
      cat("parallel nsteps=",nrow(statsmatrix),"by",
          MCMCparams.parallel$nsteps,"\n")
      ergm.rpvm.clean(rpvmbasename=rpvmbasename)
    }
  }
#   cat(paste("z$diffnwheads = ",maxchanges,z$diffnwheads[1],"\n"))
#   cat(paste("z$dissnwheads = ",maxchanges,z$dissnwheads[1],"\n"))
#   cat(paste("z$newwhead = ",maxchanges,z$newnwheads[1],"\n"))
  newnetwork<-newnw.extract(nw,z)
#   Next create the network of differences from the origianl one
  diffedgelist<-if(z$diffnwheads[1]>0){
    cbind(z$diffnwtime[2:(z$diffnwtime[1]+1)],z$diffnwheads[2:(z$diffnwheads[1]+1)],z$diffnwtails[2:(z$diffnwheads[1]+1)])
  }else{
    matrix(0, ncol=3, nrow=0)
  }
  
  colnames(statsmatrix.form) <- model.form$coef.names
  colnames(statsmatrix.diss) <- model.diss$coef.names
  
#
# recenter statsmatrix by mean statistics if necessary
#
# ms <- MCMCparams$meanstats
# print(statsmatrix[nrow(statsmatrix),])
# print(apply(statsmatrix,2,mean))
# print(ms)
# if(!is.null(ms)) {
#   if (is.null(names(ms)) && length(ms) == length(model.form$coef.names))
#     names(ms) <- model.form$coef.names
#   obs <- MCMCparams$orig.obs
#   obs <- obs[match(colnames(statsmatrix), names(obs))]
# print("adjusting:\n")
# print(obs)
#   ms  <-  ms[match(names(obs), names(ms))]
#   matchcols <- match(names(ms), names(obs))
#   if (any(!is.na(matchcols))) {
#     ms[!is.na(matchcols)] <- ms[!is.na(matchcols)] - obs[matchcols[!is.na(matchcols)]]
#     statsmatrix[,!is.na(matchcols)] <- sweep(as.matrix(
#        statsmatrix[,!is.na(matchcols)]), 2, ms[!is.na(matchcols)], "-")
#   }
# }
# print(statsmatrix[nrow(statsmatrix),])
# print(apply(statsmatrix,2,mean))

  list(statsmatrix.form=statsmatrix.form, statsmatrix.diss=statsmatrix.diss,
       newnetwork=newnetwork,
       meanstats.form=MCMCparams$meanstats.form,
       meanstats.diss=MCMCparams$meanstats.diss,
       changed=diffedgelist,
       maxchanges=MCMCparams$maxchanges)
}
