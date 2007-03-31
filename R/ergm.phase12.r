ergm.phase12 <- function(g, model,
                        MHproposal, eta0,
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
  maxedges <- max(MCMCparams$maxedges, Clist$nedges)/5
  MCMCparams$maxedges <- MCMCparams$maxedges/5
  z <- list(newnw=maxedges+1)
  while(z$newnw[1] >= maxedges){
    maxedges <- 5*maxedges
    MCMCparams$maxedges <- 5*MCMCparams$maxedges
    if(verbose){cat(paste("MCMC workspace is",maxedges,"\n"))}
#
      z <- .C("MCMCPhase12",
          as.double(Clist$heads), as.double(Clist$tails), 
          as.double(Clist$nedges), as.double(Clist$n),
          as.integer(Clist$dir), as.double(Clist$bipartite),
          as.integer(Clist$nterms), 
          as.character(Clist$fnamestring),
          as.character(Clist$snamestring),
          as.character(MHproposal$type), as.character(MHproposal$package),
          as.double(Clist$inputs),
          eta=as.double(eta0),
          as.double(MCMCparams$samplesize),
          as.double(MCMCparams$gain), as.double(ms),
          as.integer(MCMCparams$phase1),
          as.integer(MCMCparams$nsub),
          s = double(MCMCparams$samplesize * Clist$nparam),
          as.double(MCMCparams$burnin), as.double(MCMCparams$interval),
          newnw = integer(maxedges),
          as.integer(verbose), 
          as.integer(BD$attribs), 
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
    eta <- z$eta
    names(eta) <- names(eta0)
    if(z$newnw[1]>1){
      newedgelist <- matrix(z$newnw[2:z$newnw[1]], ncol=2, byrow=TRUE)
    }else{
      newedgelist <- matrix(0, ncol=2, nrow=0)
    }
    colnames(statsmatrix) <- model$coef.names
    list(statsmatrix=statsmatrix, newedgelist=newedgelist, meanstats=ms,
       maxedges=MCMCparams$maxedges,
       eta=eta)
}
