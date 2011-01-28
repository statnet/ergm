########################################################################################
# The <ergm.getMCMCsample> function samples networks using an MCMC algorithm via
# <MCMC_wrapper.C> and returns the stats matrix of the sampled networks and a single
# network as an edgelist. Note that the stats will be relative to the original network,
# i.e., the calling function must shift the statistics if required. The calling function
# must also attach column names to the statistics matrix if required.
#
# --PARAMETERS--
#   Clist     :  a list of parameters required by <MCMC_wrapper.C> and the result of
#                calling <ergm.Cprepare>
#   MHproposal:  a list of the parameters needed for Metropolis-Hastings proposals and
#                the result of calling <MHproposal>
#   eta0      :  the initial eta coefficients 
#   verbose   :  whether the C functions should be verbose; default=FALSE 
#   MCMCparams:  list of MCMC tuning parameters; those recognized include
#         maxedges      :  the maximum number of new edges that memory will be
#                          allocated for
#         samplesize    :  the number of networks to be sampled
#         nmatrixentries:  the number of entries the the returned 'statsmatrix'
#                          will have
#         interval      :  the number of proposals to ignore between sampled networks
#         burnin        :  the number of proposals to initially ignore for the burn-in
#                          period
#         Clist.miss    :  a corresponding 'Clist' for the network of missing edges,
#                          as returned by <ergm.design>
#         Clist.dt      :  yet another Clist, this one for fitting dynamic models 
#
#
# --RETURNED--
#   the sample as a list containing:
#     statsmatrix:  the stats matrix for the sampled networks, RELATIVE TO THE ORIGINAL
#                   NETWORK!
#     edgelist   :  the edgelist of the final?? sampled network
#
#########################################################################################

ergm.getMCMCsample <- function(Clist, MHproposal, eta0, MCMCparams, verbose=FALSE) {
  maxedges <- MCMCparams$maxedges
  nedges <- c(Clist$nedges,0,0)
  tails <- Clist$tails
  heads <- Clist$heads
  if(!is.null(MCMCparams$Clist.miss)){
    nedges[2] <- MCMCparams$Clist.miss$nedges
    tails <- c(tails, MCMCparams$Clist.miss$tails)
    heads <- c(heads, MCMCparams$Clist.miss$heads)
  }
  if(!is.null(MCMCparams$Clist.dt)){
    nedges[3] <- MCMCparams$Clist.dt$nedges
    tails <- c(tails, MCMCparams$Clist.dt$tails)
    heads <- c(heads, MCMCparams$Clist.dt$heads)
  }
  # *** don't forget, tails is now passed in before heads.
  z <- .C("MCMC_wrapper",
  as.integer(length(nedges)), as.integer(nedges),
  as.integer(tails), as.integer(heads),
  as.integer(Clist$maxpossibleedges), as.integer(Clist$n),
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
  newnwtails = integer(MCMCparams$maxedges),
  newnwheads = integer(MCMCparams$maxedges),
  as.integer(verbose), as.integer(MHproposal$bd$attribs),
  as.integer(MHproposal$bd$maxout), as.integer(MHproposal$bd$maxin),
  as.integer(MHproposal$bd$minout), as.integer(MHproposal$bd$minin),
  as.integer(MHproposal$bd$condAllDegExact), as.integer(length(MHproposal$bd$attribs)),
  as.integer(maxedges),
  PACKAGE="ergm")

  nedges <- z$newnwtails[1]  # This tells how many new edges there are
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
    newedgelist <- matrix(0, ncol=2, nrow=0)
  } else { 
    ## Post-processing of z$newnwtails and z$newnwheads: Combine into newedgelist
    ## The tails are listed starting at z$newnwtails[2], and similarly for heads.
    newedgelist <- cbind(z$newnwtails[2:(nedges+1)], z$newnwheads[2:(nedges+1)])
  }

  ## Post-processing of z$statsmatrix element: coerce to correct-sized matrix
  statsmatrix <- matrix(z$statsmatrix, nrow = MCMCparams$samplesize, byrow=TRUE)

  ## outta here:  Return list with "statsmatrix" and "newedgelist"
  return(list(statsmatrix = statsmatrix, newedgelist = newedgelist))
}

