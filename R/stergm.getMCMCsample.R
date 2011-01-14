############################################################################
# The <stergm.getMCMCsample> function collects a sample of networks and
# returns the formation and dissolution statistics of each sample, along with
# a toggle matrix of the changes needed from the original network to each
# in the sample
#
# --PARAMETERS--
#   nw             : a network object
#   model.form     : a formation model, as returned by <ergm.getmodel>
#   model.diss     : a dissolution model, as returned by <ergm.getmodel>
#   MHproposal.form: a list of parameters needed for MHproposals of the
#                    formations
#   MHproposal.diss: a list of parameters needed for MHproposals of the
#                    dissolutions
#   theta0.form    : the initial model coefficients for 'model.form'
#   theta0.diss    : the initial model coefficients for 'model.diss'
#   MCMCparams     : the list of parameters controlling the MCMC algorithm;
#                    recognized components include:
#      maxchanges    :  this ends up defining the "MCMCDyn workspace", but
#                       in a peculiar way: the 'maxchanges' you give, let's
#                       call this number 'mc' is fed into this equation:
#                           mc < 5^x * mc - 55
#                       the 'x' that solves the equation is then used and
#                       and 'maxchanges' and the MCMCDyn workspace is defined as
#                           5^(x-1) * (the original mc)
#                       if the computed Clist.form$nedges is > than the original
#                       'maxchanges', 'maxchanges' is ignored and this odd little
#                       process occurs to the 'nedeges'
#      meanstats.form:  the mean value parameters for 'model.form'
#      meanstats.diss:  the mean value parameters for 'model.diss'
#      parallel      :  the number of threads in which to run sampling
#      samplesize    :  the desired MCMC sample size
#      toggles       :  whether to return the toggle matrix; non-NULL values
#                       will include 'changed' in the return list, NULL
#                       will not
#      MH.burnin     :  this is accepted as 'MH_interval' which determines
#                       the number of proposals in each MCMC step
#      time.burnin   :  the number of MCMC steps to disregard for the burnin
#                       period 
#      time.interval :  the number of MCMC steps to disregard in between
#                       sampled networks
#   verbose        : whether this and other functions should be verbose
#
# Note:  In reality, there should be many fewer arguments to this function,
# since most info should be passed via Clist (this is, after all, what Clist
# is for:  Holding all arguments required for the .C call).  In particular,
# the elements of MHproposal.form, MCMCparams, verbose, should certainly
# be part of Clist.  But this is a project for another day!
#
# --RETURNED--
#   the MCMC sample as a list containing:
#     statsmatrix.form: the matrix of sampled statistics for 'model.form'
#     statsmatrix.diss: the matrix of sampled statistics for 'model.form'
#     newnetwork      : the final network from the sampling process
#     meanstats.form  : the 'meanstats.form' passed in via 'MCMCparams' 
#     meanstats.diss  : the 'meanstats.diss' passed in via 'MCMCparams' 
#     changed         : a toggle matrix, where the first column is
#                       the timestamp of the toggle and the 2nd and 3rd
#                       columns are the head & tail of the toggle; this
#                       is only returned if 'MCMCparams'$toggles is not NULL
#     maxchanges      : the "MCMC Dyn workspace"; see 'maxchanges' in the
#                       input param list
#
############################################################################
  
stergm.getMCMCsample <- function(nw, model.form, model.diss,
                                  MHproposal.form, MHproposal.diss, theta.form0, theta.diss0, MCMCparams, 
                                  verbose){
#
#   Check for truncation of the returned edge list
#
  Clist.form <- ergm.Cprepare(nw, model.form)
  Clist.diss <- ergm.Cprepare(nw, model.diss)

  ## Sanity check: parameter vectors are as long as they are supposed to be:

  if(length(theta.form0)!=model.form$etamap$etalength) stop("Wrong number of formation parameters for the model: ",length(theta.form0), ", but should be ", model.form$etamap$etalength)
  if(length(theta.diss0)!=model.diss$etamap$etalength) stop("Wrong number of dissolution parameters for the model: ",length(theta.diss0), ", but should be ", model.diss$etamap$etalength)
  
  maxchanges <- max(MCMCparams$maxchanges, Clist.form$nedges)/5
  MCMCparams$maxchanges <- MCMCparams$maxchanges/5
  if(is.null(MCMCparams$meanstats.form)) MCMCparams$meanstats.form<-numeric(length(model.form$coef.names))
  if(is.null(MCMCparams$meanstats.diss)) MCMCparams$meanstats.diss<-numeric(length(model.diss$coef.names))
  z <- list(newnwheads=maxchanges+1,diffnwheads=maxchanges+1)
  while(z$newnwheads[1]  >= maxchanges - 10 || 
        z$diffnwheads[1] >= maxchanges - 10){
    maxchanges <- 5*maxchanges
    MCMCparams$maxchanges <- 5*MCMCparams$maxchanges
    if(verbose){cat(paste("MCMCDyn workspace is",maxchanges,"\n"))}
    
    if(MCMCparams$parallel==0){
      z <- .C("MCMCDyn_wrapper",
              # Observed network.
              as.integer(Clist.form$heads), as.integer(Clist.form$tails), 
              as.integer(Clist.form$nedges), as.integer(Clist.form$maxpossibleedges),
              as.integer(Clist.form$n),
              as.integer(Clist.form$dir), as.integer(Clist.form$bipartite),
              # Ordering of formation and dissolution.
              as.integer(Clist.diss$stergm.order.code),
              # Formation terms and proposals.
              as.integer(Clist.form$nterms), 
              as.character(Clist.form$fnamestring),
              as.character(Clist.form$snamestring),
              as.character(MHproposal.form$name), as.character(MHproposal.form$package),
              as.double(Clist.form$inputs), as.double(theta.form0),
              #?  Add:  as.double(length(MHproposal.form$args)), as.double(MHproposal.form$args),
              # Dissolution terms and proposals.
              as.integer(Clist.diss$nterms), 
              as.character(Clist.diss$fnamestring),
              as.character(Clist.diss$snamestring),
              as.character(MHproposal.diss$name), as.character(MHproposal.diss$package),
              as.double(Clist.diss$inputs), as.double(theta.diss0),
              # Degree bounds.
              as.integer(MHproposal.form$bd$attribs), 
              as.integer(MHproposal.form$bd$maxout), as.integer(MHproposal.form$bd$maxin),
              as.integer(MHproposal.form$bd$minout), as.integer(MHproposal.form$bd$minin),
              as.integer(MHproposal.form$bd$condAllDegExact), as.integer(length(MHproposal.form$bd$attribs)),
              # MCMC settings.
              as.double(MCMCparams$samplesize), as.integer(MCMCparams$MH.burnin),
              as.double(MCMCparams$time.burnin), as.double(MCMCparams$time.interval),
              # Space for output.
              s.form = as.double(cbind(MCMCparams$meanstats.form,matrix(0,nrow=length(model.form$coef.names),ncol=MCMCparams$samplesize))),
              s.diss = as.double(cbind(MCMCparams$meanstats.diss,matrix(0,nrow=length(model.diss$coef.names),ncol=MCMCparams$samplesize))),
              newnwheads = integer(maxchanges), newnwtails = integer(maxchanges), 
              as.double(maxchanges),
              diffnwtime = integer(maxchanges),
              diffnwheads = integer(maxchanges),
              diffnwtails = integer(maxchanges),
              as.integer(verbose), 
              PACKAGE="ergm") 
      statsmatrix.form <- matrix(z$s.form, nrow=MCMCparams$samplesize+1,
                                 ncol=Clist.form$nstats,
                                 byrow = TRUE)[-1,,drop=FALSE]
      statsmatrix.diss <- matrix(z$s.diss, nrow=MCMCparams$samplesize+1,
                                 ncol=Clist.diss$nstats,
                                 byrow = TRUE)[-1,,drop=FALSE]
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
      outsetuppvm <- stergm.rpvm.setup(rpvmbasename, verbose=verbose,
                                       packagename=packagename)
#
#   Saving the common variables
#
      save(Clist.form,
           Clist.diss,
           MHproposal.form,
           MHproposal.diss,
           theta.form0,
           theta.diss0,
           MCMCparams.parallel,
           maxchanges, 
           verbose,
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
                                  matrix(z$s, nrow=MCMCparams.parallel$samplesize,
                                         ncol=Clist.form$nstats,
                                         byrow = TRUE))
        statsmatrix.diss <- rbind(statsmatrix.diss,
                                      matrix(z$s, nrow=MCMCparams.parallel$samplesize,
                                             ncol=Clist.diss$nstats,
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
#   cat(paste("z$diffnwheads = ",maxchanges,z$diffnwheads[1],"\n"))
#   cat(paste("z$dissnwheads = ",maxchanges,z$dissnwheads[1],"\n"))
#   cat(paste("z$newwhead = ",maxchanges,z$newnwheads[1],"\n"))
  newnetwork<-newnw.extract(nw,z)
#   Next create the network of differences from the origianl one

  
  diffedgelist<-if(!is.null(MCMCparams$toggles)) {
    if(z$diffnwheads[1]>0){
      cbind(z$diffnwtime[2:(z$diffnwtime[1]+1)],z$diffnwheads[2:(z$diffnwheads[1]+1)],z$diffnwtails[2:(z$diffnwheads[1]+1)])
    }else{
      matrix(0, ncol=3, nrow=0)
    }
  }else{
    NULL
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
