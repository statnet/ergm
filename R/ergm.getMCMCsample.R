# This function is nothing other than an R wrapper for the MCMC_wrapper
# function in C.  It assumes that the calling function will send only what is 
# necessary, namely:
#
#  Clist (the result of calling ergm.Cprepare(network, model))
#  MHproposal (the result of calling MHproposal(...))
#  eta0 (the canonical parameter value governing the MCMC simulation)
#  MCMCparams (sort of a catch-all for other arguments passed)
#  verbose (which governs the verbosity of the C functions)
#
#  It returns only the named elements of the .C() call, after some 
#  post-processing (e.g., the statistics matrix is coerced to the correct
#                   dimensions and given appropriate column names)
#  NB:  The statistics are all RELATIVE TO THE ORIGINAL MATRIX!
#       i.e., the calling function must shift the statistics if necessary.

ergm.getMCMCsample <- function(Clist, MHproposal, eta0, MCMCparams, verbose=FALSE) {

  if(verbose)
    cat("\nheads: ",
    as.integer(Clist$heads), 
    "\ntails: ",
    as.integer(Clist$tails),
    "\nnedges: ",
    as.integer(Clist$nedges), 
    "\nmaxpossibleedges: ",
    as.integer(Clist$maxpossibleedges), 
    "\nn: ",
    as.integer(Clist$n),
    "\ndir: ",
    as.integer(Clist$dir), 
    "\nbipartite: ",
    as.integer(Clist$bipartite),
    "\nnterms: ",
    as.integer(Clist$nterms),
    "\nfnamestring: ",
    as.character(Clist$fnamestring),
    "\nsnamestring: ",
    as.character(Clist$snamestring),
    "\nMHproposalname: ",
    as.character(MHproposal$name), 
    "\nMHproposalpackage: ",
    as.character(MHproposal$package),
    "\ninputs: ",
    as.double(Clist$inputs), 
    "\neta0: ",
    as.double(eta0),
    "\nsamplesize: ",
    as.integer(MCMCparams$samplesize),
    "\nstats: ",
    #  s = as.double(t(MCMCparams$stats)),
    double(MCMCparams$nmatrixentries),
    "\nburnin: ",
    as.integer(MCMCparams$burnin), 
    "\ninterval: ",
    as.integer(MCMCparams$interval),
    "\nnumnewheads and numnewtails same as maxedges",
    #  newnwheads = integer(maxedges),
    #  newnwtails = integer(maxedges),
    "\nverbose: ",
    as.integer(verbose), 
    "\nbd attribs: ",
    as.integer(MHproposal$bd$attribs),
    "\nbd maxout: ",
    as.integer(MHproposal$bd$maxout), 
    "\nbd maxin: ",
    as.integer(MHproposal$bd$maxin),
    "\nbd minout: ",
    as.integer(MHproposal$bd$minout), 
    "\nbf condAllDegExact: ",
    as.integer(MHproposal$bd$minin),
    "\nbd condAllDegExact: ",
    as.integer(MHproposal$bd$condAllDegExact), 
    "\nbd attribs: ",
    as.integer(length(MHproposal$bd$attribs)),
    "\nmaxedges: ",
    as.integer(MCMCparams$maxedges),
    "\nmiss.heads: ",
    as.integer(MCMCparams$Clist.miss$heads), 
    "\nmiss.tails: ",
    as.integer(MCMCparams$Clist.miss$tails),
    "\nmiss.nedges: ",
    as.integer(MCMCparams$Clist.miss$nedges))


  z <- .C("MCMC_wrapper",
  as.integer(Clist$heads), as.integer(Clist$tails),
  as.integer(Clist$nedges), as.integer(Clist$maxpossibleedges), as.integer(Clist$n),
  as.integer(Clist$dir), as.integer(Clist$bipartite),
  as.integer(Clist$nterms),
  as.character(Clist$fnamestring),
  as.character(Clist$snamestring),
  as.character(MHproposal$name), as.character(MHproposal$package),
  as.double(Clist$inputs), as.double(eta0),
  as.integer(MCMCparams$samplesize),
  # The line below was changed as of version 2.2-3.  Now, the statsmatrix is 
  # initialized to zero instead of allowing the first row to be nonzero, then 
  # adding this first row to each row within MCMC_wrapper.
  # Any unmodified old function trying to use the new version will generate an 
  # error because the MCMCparams$nmatrixentries object is new and will not yet 
  # exist in an unmodified function.
  statsmatrix = double(MCMCparams$nmatrixentries),
  #  statsmatrix = as.double(t(MCMCparams$stats)), # By default, as.double goes bycol, not byrow; thus, we use the transpose here.
  as.integer(MCMCparams$burnin), 
  as.integer(MCMCparams$interval),
  newnwheads = integer(MCMCparams$maxedges),
  newnwtails = integer(MCMCparams$maxedges),
  as.integer(verbose), as.integer(MHproposal$bd$attribs),
  as.integer(MHproposal$bd$maxout), as.integer(MHproposal$bd$maxin),
  as.integer(MHproposal$bd$minout), as.integer(MHproposal$bd$minin),
  as.integer(MHproposal$bd$condAllDegExact), as.integer(length(MHproposal$bd$attribs)),
  as.integer(MCMCparams$maxedges),
  as.integer(MCMCparams$Clist.miss$heads), as.integer(MCMCparams$Clist.miss$tails),
  as.integer(MCMCparams$Clist.miss$nedges),
  PACKAGE="ergm")

  ## Post-processing of z$statsmatrix element: coerce to correct-sized matrix
  statsmatrix <- matrix(z$statsmatrix, nrow = MCMCparams$samplesize, byrow=TRUE)
  
  ## Post-processing of z$newnwheads and z$newnwtails: Combine into newedgelist
  nedges <- z$newnwheads[1]  # This tells how many new edges there are, whose
       # heads are listed starting at z$newnwheads[2], and similarly for tails.
  if (nedges==0) { newedgelist <- matrix(0, ncol=2, nrow=0)}
  else { newedgelist <- cbind(z$newnwtails[2:(nedges+1)], z$newnwheads[2:(nedges+1)])}
  
  ## outta here:  Return list with "statsmatrix" and "newedgelist"
  return(list(statsmatrix = statsmatrix, newedgelist = newedgelist))
}


# ergm.mcmcslave should now be unnecessary; the ergm.getMCMCsample function
# now takes its place (as of version 2.2-3)

# Function the slaves will call to perform a validation on the
# mcmc equal to their slave number.
# Assumes: Clist MHproposal eta0 MCMCparams maxedges verbose
ergm.mcmcslave <- function(Clist,MHproposal,eta0,MCMCparams,maxedges,verbose) {
  warning("Using deprecated function ergm.mcmcslave!")
  z <- .C("MCMC_wrapper",
  as.integer(Clist$heads), as.integer(Clist$tails),
  as.integer(Clist$nedges), as.integer(Clist$maxpossibleedges), as.integer(Clist$n),
  as.integer(Clist$dir), as.integer(Clist$bipartite),
  as.integer(Clist$nterms),
  as.character(Clist$fnamestring),
  as.character(Clist$snamestring),
  as.character(MHproposal$name), as.character(MHproposal$package),
  as.double(Clist$inputs), as.double(eta0),
  as.integer(MCMCparams$samplesize),
#  s = as.double(t(MCMCparams$stats)),
  s = double(MCMCparams$samplesize * length(MCMCparams$stats)),
  as.integer(MCMCparams$burnin), 
  as.integer(MCMCparams$interval),
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
  # save the results
  list(s=z$s, newnwheads=z$newnwheads, newnwtails=z$newnwtails)
}
