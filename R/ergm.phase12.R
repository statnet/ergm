#  File ergm/R/ergm.phase12.R
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
ergm.phase12 <- function(g, model,
                        MHproposal, eta0,
                        MCMCparams, verbose) {
# ms <- MCMCparams$meanstats
# if(!is.null(ms)) {
#   if (is.null(names(ms)) && length(ms) == length(model$coef.names))
#     names(ms) <- model$coef.names
#   obs <- MCMCparams$orig.obs
#   obs <- obs[match(names(ms), names(obs))]
#   ms  <-  ms[match(names(obs), names(ms))]
#   matchcols <- match(names(ms), names(obs))
#   if (any(!is.na(matchcols))) {
#     ms[!is.na(matchcols)] <- ms[!is.na(matchcols)] - obs[matchcols[!is.na(matchcols)]]
#   }
# }
  Clist <- ergm.Cprepare(g, model)
  maxedges <- max(MCMCparams$maxedges, Clist$nedges)/5
  MCMCparams$maxedges <- MCMCparams$maxedges/5
  z <- list(newnwheads=maxedges+1)
  while(z$newnwheads[1] >= maxedges){
    maxedges <- 5*maxedges
    MCMCparams$maxedges <- 5*MCMCparams$maxedges
    if(verbose){cat(paste("MCMC workspace is",maxedges,"\n"))}
#
    z <- .C("MCMCPhase12",
            as.integer(Clist$heads), as.integer(Clist$tails), 
            as.integer(Clist$nedges), as.integer(Clist$maxpossibleedges),
            as.integer(Clist$n),
            as.integer(Clist$dir), as.integer(Clist$bipartite),
            as.integer(Clist$nterms), 
            as.character(Clist$fnamestring),
            as.character(Clist$snamestring),
            as.character(MHproposal$name), as.character(MHproposal$package),
            as.double(Clist$inputs),
            eta=as.double(eta0),
            as.integer(MCMCparams$samplesize),
            as.double(MCMCparams$gain), as.double(MCMCparams$stats),
            as.integer(MCMCparams$phase1),
            as.integer(MCMCparams$nsub),
            s = double(MCMCparams$samplesize * Clist$nstats),
            as.integer(MCMCparams$burnin), as.integer(MCMCparams$interval),
            newnwheads = integer(maxedges),
            newnwtails = integer(maxedges),
            as.integer(verbose), 
            as.integer(MHproposal$bd$attribs), 
            as.integer(MHproposal$bd$maxout), as.integer(MHproposal$bd$maxin),
            as.integer(MHproposal$bd$minout), as.integer(MHproposal$bd$minin),
            as.integer(MHproposal$bd$condAllDegExact), as.integer(length(MHproposal$bd$attribs)), 
            as.integer(maxedges),
            as.integer(0.0), as.integer(0.0), 
            as.integer(0),
            PACKAGE="ergm") 
  }
  statsmatrix <- matrix(z$s, nrow=MCMCparams$samplesize,
                        ncol=Clist$nstats,
                        byrow = TRUE)
  eta <- z$eta
  names(eta) <- names(eta0)

  newnetwork<-newnw.extract(g,z)
  
  colnames(statsmatrix) <- model$coef.names
  list(statsmatrix=statsmatrix, newnetwork=newnetwork, meanstats=MCMCparams$meanstats,
       maxedges=MCMCparams$maxedges,
       eta=eta)
}
