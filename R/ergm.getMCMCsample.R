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
#  post-processing: the statistics matrix is coerced to the correct
#                   dimensions and the heads/tails are returned as an edgelist
#  NB:  The statistics are all RELATIVE TO THE ORIGINAL NETWORK!
#       i.e., the calling function must shift the statistics if required.
#       The calling function must also attach column names to the statistics
#       matrix if required.

ergm.getMCMCsample <- function(Clist, MHproposal, eta0, MCMCparams, verbose=FALSE) {
  maxedges <- MCMCparams$maxedges
  if(is.null(Clist$weights))
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
  # exist in an unmodified function.  This is worth it:  There is no reason
  # that MCMCparams should include a huge matrix.
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
            as.integer(maxedges),
            as.integer(MCMCparams$Clist.miss$heads), as.integer(MCMCparams$Clist.miss$tails),
            as.integer(MCMCparams$Clist.miss$nedges),
            PACKAGE="ergm")
  else
    z <- .C("WtMCMC_wrapper",
            as.integer(Clist$heads), as.integer(Clist$tails), as.double(Clist$weights),
            as.integer(Clist$nedges), as.double(Clist$baseline_weight), as.integer(Clist$maxpossibleedges), as.integer(Clist$n),
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
  # exist in an unmodified function.  This is worth it:  There is no reason
  # that MCMCparams should include a huge matrix.
            statsmatrix = double(MCMCparams$nmatrixentries),
  #  statsmatrix = as.double(t(MCMCparams$stats)), # By default, as.double goes bycol, not byrow; thus, we use the transpose here.
            as.integer(MCMCparams$burnin),
            as.integer(MCMCparams$interval),
            newnwheads = integer(MCMCparams$maxedges),
            newnwtails = integer(MCMCparams$maxedges),
            newnwweights = double(MCMCparams$maxedges),
            as.integer(verbose), 
            as.integer(maxedges),
            PACKAGE="ergm")


  nedges <- z$newnwheads[1]  # This tells how many new edges there are
  if (nedges >= maxedges) {
    # The simulation has filled up the available memory for storing edges, 
    # so rerun it with ten times more
    # To do:  Check to see whether it is possible to pass a "statsonly"
    # argument to the C code, thus avoiding the need to store the final network
    # and eliminating the need to make this particular check.
    MCMCparams$maxedges <- maxedges * 10
    if (verbose) cat("Increasing possible number of newedges to ", 
                     MCMCparams$maxedges, "\n")
    return(ergm.getMCMCsample(Clist, MHproposal, eta0, MCMCparams, verbose=FALSE))
  } else if (nedges==0) { 
    newedgelist <- matrix(0, ncol=2+(!is.null(Clist$weights)), nrow=0)
  } else { 
    ## Post-processing of z$newnwheads and z$newnwtails: Combine into newedgelist
    ## The heads are listed starting at z$newnwheads[2], and similarly for tails.
    newedgelist <- cbind(z$newnwtails[2:(nedges+1)], z$newnwheads[2:(nedges+1)])
    if(!is.null(Clist$weights)) newedgelist<-cbind(newedgelist,z$newnwweights[2:(nedges+1)])
  }

  ## Post-processing of z$statsmatrix element: coerce to correct-sized matrix
  statsmatrix <- matrix(z$statsmatrix, nrow = MCMCparams$samplesize, byrow=TRUE)

  ## outta here:  Return list with "statsmatrix" and "newedgelist"
  return(list(statsmatrix = statsmatrix, newedgelist = newedgelist))
}

