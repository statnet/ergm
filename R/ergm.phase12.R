#  File R/ergm.phase12.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
###############################################################################
# The <ergm.phase12> function is a wrapper for the <MCMC.phase12.C> method,
# which collects a sample of networks and returns the matrix of summary
# statistics
#
# --PARAMETERS--
#   g         : a network object
#   model     : a model for 'g', as returned by <ergm.getmodel>
#   MHproposal: an MHproposal object, as returned by <MHproposal>
#   eta0      : the vector of initial eta coefficients
#   control: a list of control parameters for the MCMC algorithm;
#               recognized components include:
#                     'maxedges'     'samplesize'     'gain'
#                     'stats'        'phase1'         'nsub'
#                     'burnin'       'interval'       'target.stats'
#               the purpose of most of these variables is given in the
#               <control.ergm> function header; 'stats' seems to be
#                used as the mean statistics; 'target.stats' is merely
#                returned.
#   verbose   : whether the C functions should be verbose (T or F)
#
# --RETURNED--
#   a list containing
#     statsmatrix: the matrix of summary statistics
#     newnetwork : the final network sampled
#     target.stats  : the 'target.stats' from 'control'
#     maxedges   : the 'maxedges' from 'control'
#     eta        : the parameters used to produce the sample given
#                  by 'statsmatrix'
#
###############################################################################

ergm.phase12 <- function(g, model,
                        MHproposal, eta0,
                        control, verbose) {
# ms <- model$target.stats
# if(!is.null(ms)) {
#   if (is.null(names(ms)) && length(ms) == length(model$coef.names))
#     names(ms) <- model$coef.names
#   obs <- control$orig.obs
#   obs <- obs[match(names(ms), names(obs))]
#   ms  <-  ms[match(names(obs), names(ms))]
#   matchcols <- match(names(ms), names(obs))
#   if (any(!is.na(matchcols))) {
#     ms[!is.na(matchcols)] <- ms[!is.na(matchcols)] - obs[matchcols[!is.na(matchcols)]]
#   }
# }
  Clist <- ergm.Cprepare(g, model)
  maxedges <- max(control$MCMC.init.maxedges, Clist$nedges)/5
  control$MCMC.init.maxedges <- control$MCMC.init.maxedges/5
  z <- list(newnwtails=maxedges+1)
  while(z$newnwtails[1] >= maxedges){
    maxedges <- 5*maxedges
    control$MCMC.init.maxedges <- 5*control$MCMC.init.maxedges
    if(verbose){cat(paste("MCMC workspace is",maxedges,"\n"))}
    # *** don't forget, pass in tails first now, not heads
    z <- .C("MCMCPhase12",
            as.integer(Clist$tails), as.integer(Clist$heads), 
            as.integer(Clist$nedges), 
            as.integer(Clist$n),
            as.integer(Clist$dir), as.integer(Clist$bipartite),
            as.integer(Clist$nterms), 
            as.character(Clist$fnamestring),
            as.character(Clist$snamestring),
            as.character(MHproposal$name), as.character(MHproposal$pkgname),
            as.double(Clist$inputs),
            eta=as.double(eta0),
            as.integer(control$MCMC.samplesize),
            as.double(control$gain), as.double(control$stats),
            as.integer(control$phase1),
            as.integer(control$nsub),
            s = double(control$MCMC.samplesize * Clist$nstats),
            as.integer(control$MCMC.burnin), as.integer(control$MCMC.interval),
            newnwtails = integer(maxedges),
            newnwheads = integer(maxedges),
            as.integer(verbose), 
            as.integer(MHproposal$arguments$constraints$bd$attribs), 
            as.integer(MHproposal$arguments$constraints$bd$maxout), as.integer(MHproposal$arguments$constraints$bd$maxin),
            as.integer(MHproposal$arguments$constraints$bd$minout), as.integer(MHproposal$arguments$constraints$bd$minin),
            as.integer(MHproposal$arguments$constraints$bd$condAllDegExact), as.integer(length(MHproposal$arguments$constraints$bd$attribs)), 
            as.integer(maxedges),
            as.integer(0.0), as.integer(0.0), 
            as.integer(0),
            PACKAGE="ergm") 
  }
  statsmatrix <- matrix(z$s, nrow=control$MCMC.samplesize,
                        ncol=Clist$nstats,
                        byrow = TRUE)
   eta <- z$eta
  names(eta) <- names(eta0)

  newnetwork<-newnw.extract(g,z)
  
  colnames(statsmatrix) <- model$coef.names
  list(statsmatrix=statsmatrix, newnetwork=newnetwork, target.stats=model$target.stats,
       maxedges=control$MCMC.init.maxedges,
       eta=eta)
}
