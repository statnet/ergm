#  File ergm/R/ergm.getMCMCsample.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of Washington
#                David R. Hunter, Penn State University
#                Carter T. Butts, University of California - Irvine
#                Steven M. Goodreau, University of Washington
#                Martina Morris, University of Washington
# Copyright 2007 The statnet Development Team
######################################################################
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
#  Parallel running
#
    if(MCMCparams$parallel==0){
    flush.console()
    z <- ergm.mcmcslave(Clist,MHproposal,eta0,MCMCparams,maxedges,verbose)
    nedges <- z$newnwheads[1]
    statsmatrix <- matrix(z$s, nrow=MCMCparams$samplesize,
                          ncol=Clist$nstats,
                          byrow = TRUE)
    newnetwork <- newnw.extract(nw,z)
    if(nedges >= 50000-1){
      cat("\n Warning:")
      cat("\n   The network has more than 50000 edges, and the model is likely to be degenerate.\n")
#  NOT SURE ABOUT COMMENTING OUT THE FOLLOWING THREE LINES:
#      statsmatrix <- matrix(0, nrow=MCMCparams$samplesize,
#                            ncol=Clist$nstats)
#      newnetwork <- nw
    }      
    }else{
      stop("parallization not enabled for now.")
    }
  }
  colnames(statsmatrix) <- model$coef.names

##
## recenter statsmatrix by mean statistics if necessary
##
#   ms <- Clist$meanstats
#   if(!is.null(ms)) {
#     if (is.null(names(ms)) && length(ms) == length(model$coef.names))
#       names(ms) <- model$coef.names
##    obs <- summary(model$formula)
#     obs <- Clist$obs
##    print(paste("obs=",obs))
##    print(paste("statsmatrix=",apply(statsmatrix,2,mean)))
#     obs <- obs[match(colnames(statsmatrix), names(obs))]
#     ms  <-  ms[match(names(obs), names(ms))]
#     matchcols <- match(names(ms), names(obs))
#     if (any(!is.na(matchcols))) {
#       ms[!is.na(matchcols)] <- ms[!is.na(matchcols)] - obs[matchcols[!is.na(matchcols)]]
#       statsmatrix[,!is.na(matchcols)] <- sweep(as.matrix(
#          statsmatrix[,!is.na(matchcols)]), 2, ms[!is.na(matchcols)], "-")
#     }
#   }
  list(statsmatrix=statsmatrix, newnetwork=newnetwork, 
       meanstats=Clist$meanstats, nedges=nedges)
}
# Function the slaves will call to perform a validation on the
# mcmc equal to their slave number.
# Assumes: Clist MHproposal eta0 MCMCparams maxedges verbose
ergm.mcmcslave <- function(Clist,MHproposal,eta0,MCMCparams,maxedges,verbose) {
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
  as.integer(MCMCparams$Clist.miss$heads), as.integer(MCMCparams$Clist.miss$tails),
  as.integer(MCMCparams$Clist.miss$nedges),
  PACKAGE="ergm")
  # save the results
  list(s=z$s, newnwheads=z$newnwheads, newnwtails=z$newnwtails)
}
