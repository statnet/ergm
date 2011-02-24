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
#   MCMCparams: a list of control parameters for the MCMC algorithm;
#               recognized components include:
#                     'maxedges'     'samplesize'     'gain'
#                     'stats'        'phase1'         'nsub'
#                     'burnin'       'interval'       'meanstats'
#               the purpose of most of these variables is given in the
#               <control.ergm> function header; 'stats' seems to be
#                used as the mean statistics; 'meanstats' is merely
#                returned.
#   verbose   : whether the C functions should be verbose (T or F)
#
# --RETURNED--
#   a list containing
#     statsmatrix: the matrix of summary statistics
#     newnetwork : the final network sampled
#     meanstats  : the 'meanstats' from 'MCMCparams'
#     maxedges   : the 'maxedges' from 'MCMCparams'
#     eta        : the parameters used to produce the sample given
#                  by 'statsmatrix'
#
###############################################################################

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
  z <- list(newnwtails=maxedges+1)
  while(z$newnwtails[1] >= maxedges){
    maxedges <- 5*maxedges
    MCMCparams$maxedges <- 5*MCMCparams$maxedges
    if(verbose){cat(paste("MCMC workspace is",maxedges,"\n"))}
    # *** don't forget, pass in tails first now, not heads
    z <- .C("MCMCPhase12",
            as.integer(Clist$tails), as.integer(Clist$heads), 
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
            newnwtails = integer(maxedges),
            newnwheads = integer(maxedges),
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
