ergm.getMCMCsample <- function(nw, model, MHproposal, eta0, MCMCparams, 
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
    flush.console()
    z <- ergm.mcmcslave(Clist,MHproposal,eta0,MCMCparams,maxedges,verbose)
    nedges <- z$newnwheads[1]
    if(nedges >= 50000-1){
      cat("\n Warning:")
      cat("\n   The network has more than 50000 edges, and the model is likely to be degenerate.\n")
      statsmatrix <- matrix(0, nrow=MCMCparams$samplesize,
                            ncol=Clist$nparam)
      newnetwork <- nw
    }else{
      statsmatrix <- matrix(z$s, nrow=MCMCparams$samplesize,
                            ncol=Clist$nparam,
                            byrow = TRUE)

      newnetwork <- newnw.extract(nw,z)
    }
  }
  colnames(statsmatrix) <- model$coef.names

  list(statsmatrix=statsmatrix,
       newnetwork=newnetwork,
       meanstats=Clist$meanstats,
       nedges=nedges)
}

ergm.mcmcslave <- function(Clist,MHproposal,eta0,MCMCparams,maxedges,verbose) {
  z <- .C("MCMC_wrapper",
  as.integer(Clist$heads), as.integer(Clist$tails),
  as.integer(Clist$nedges), as.integer(Clist$n),
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
  as.integer(0.0), as.integer(0.0),
  as.integer(0.0),
  PACKAGE="ergm")
  # save the results
  list(s=z$s, newnwheads=z$newnwheads, newnwtails=z$newnwtails)
}
