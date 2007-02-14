ergm.getMCMCsample <- function(Clist, model, MHproposal, eta0, MCMCparams, 
                               verbose, BD) {
# Note:  In reality, there should be many fewer arguments to this function,
# since most info should be passed via Clist (this is, after all, what Clist
# is for:  Holding all arguments required for the .C call).  In particular,
# the elements of MHproposal, MCMCparams, verbose, and BD should certainly
# be part of Clist.  But this is a project for another day!
  maxedges <- max(5000, Clist$nedges)
#
#   Check for truncation of the returned edge list
#
  z <- list(newnw=maxedges+1)
  while(z$newnw[1] > maxedges){
   maxedges <- 10*maxedges
   z <- .C("MCMC_wrapper",
          as.double(Clist$heads), as.double(Clist$tails), 
          as.double(Clist$nedges), as.double(Clist$n),
          as.integer(Clist$dir), as.double(Clist$bipartite),
          as.integer(Clist$nterms), 
          as.character(Clist$fnamestring),
          as.character(Clist$snamestring),
          as.character(MHproposal$type), as.character(MHproposal$package),
          as.double(Clist$inputs), as.double(eta0),
          as.double(MCMCparams$samplesize),
          s = double(MCMCparams$samplesize * Clist$nparam),
          as.double(MCMCparams$burnin), as.double(MCMCparams$interval), 
          newnw = integer(maxedges), 
          as.integer(verbose), as.integer(BD$attribs), 
          as.integer(BD$maxout), as.integer(BD$maxin),
          as.integer(BD$minout), as.integer(BD$minin),
          as.integer(BD$condAllDegExact), as.integer(length(BD$attribs)), 
          as.double(maxedges),
          as.double(0.0), as.double(0.0), 
          as.double(0.0), as.integer(0),
          PACKAGE="statnet") 
  }
  statsmatrix <- matrix(z$s, nrow=MCMCparams$samplesize,
                        ncol=Clist$nparam,
                        byrow = TRUE)
  colnames(statsmatrix) <- model$coef.names
  if(z$newnw[1]>1){
    newedgelist <- matrix(z$newnw[2:z$newnw[1]], ncol=2, byrow=TRUE)
  }else{
    newedgelist <- matrix(0, ncol=2, nrow=0)
  }

#
# recenter statsmatrix by mean statistics if necessary
#
    ms <- Clist$meanstats
    if(!is.null(ms)) {
      if (is.null(names(ms)) && length(ms) == length(model$coef.names))
        names(ms) <- model$coef.names
#     obs <- summary(model$formula)
      obs <- Clist$obs
#     print(paste("obs=",obs))
#     print(paste("statsmatrix=",apply(statsmatrix,2,mean)))
      obs <- obs[match(colnames(statsmatrix), names(obs))]
      ms  <-  ms[match(names(obs), names(ms))]
      matchcols <- match(names(ms), names(obs))
      if (any(!is.na(matchcols))) {
        ms[!is.na(matchcols)] <- ms[!is.na(matchcols)] - obs[matchcols[!is.na(matchcols)]]
        statsmatrix[,!is.na(matchcols)] <- sweep(as.matrix(
           statsmatrix[,!is.na(matchcols)]), 2, ms[!is.na(matchcols)], "-")
      }
    }
  list(statsmatrix=statsmatrix, newedgelist=newedgelist, meanstats=ms)
}
