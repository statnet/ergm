#===========================================================================
# This file contains the 2 following functions for ??
#       <ergm.runPILAsampler>
#       <ergm.PILAslave>
#===========================================================================




###########################################################################
# The <ergm.runPILAsampler> collects and returns a ?? sample 
#
# --PARAMETERS--
#   nw        : a network object
#   model     : a model for 'nw', as returned by <ergm.getmodel>
#   MHproposal: an MHproposal object, as returned by <MHproposal>
#   eta0      : the vector of initial eta coefficients
#   MCMCparams: a list of control parameters for the MCMC algorithm;
#               recognized components include:
#     samplesize     : the size of the sample to collect
#     interval       : the number of proposals to disregard between
#                      samples
#     stats          : ??
#     burnin         : the number of proposals to initially disregard
#                      for the burn-in period
#     PILA.steplength: ??
#     PILA.gamma     : ??
#     Clist.miss     : a list of parameters for the missing subgraph
#                      of 'nw' as returned by <ergm.design>
#   verbose   : whether the C functions should be verbose (T or F)                  
#   
# --RETURNED--
#   out: the sample as a list containing:
#     statsmatrix: the stats matrix
#     etamatrix  : the ??
#
###########################################################################

ergm.runPILAsampler <- function(nw, model, MHproposal, eta0, MCMCparams, 
                               verbose) {

  Clist <- ergm.Cprepare(nw, model)
  z <- ergm.PILAslave(Clist,MHproposal,eta0,MCMCparams,verbose,debug=verbose>0)
  nedges <- z$newnwheads[1]
  out<-list(statsmatrix=matrix(z$s, nrow=MCMCparams$samplesize,
              ncol=Clist$nstats,
              byrow = TRUE),
            etamatrix=matrix(z$eta, nrow=MCMCparams$samplesize,
              ncol=Clist$nstats,
              byrow=TRUE)
            )

  z$s<-z$eta<-NULL
  for(name in names(z)){
    out[[name]]<-matrix(z[[name]],nrow=MCMCparams$samplesize+1,byrow=TRUE)
  }

  colnames(out$statsmatrix) <- model$coef.names

  out
}




#############################################################################
# The <ergm.PILAslave> function ?? via <PILA_wrapper.C>
# Function the slaves will call to perform a validation on the
# mcmc equal to their slave number.
#
# --PARAMETERS--
#   (same as the function above, but also includes)
#   Clist: the list of items needed by the C code and returned by
#          <ergm.Cprepare>
#   debug: whether to return several of the intermediary data vectors used
#          by <PILA_wrapper.C> for debugging purposes; default=FALSE
#
# --RETURNED--
#   z: the PILA sample as a list containing:
#        s  : the stats matrix
#        eta: the ??
#      and optionally including the 9 additional debugging vectors
#
#############################################################################

ergm.PILAslave <- function(Clist,MHproposal,eta0,MCMCparams,verbose,debug=FALSE) {
  z <- .C("PILA_wrapper",
          as.integer(Clist$heads), as.integer(Clist$tails),
          as.integer(Clist$nedges),
          as.integer(Clist$n), as.integer(Clist$dir), as.integer(Clist$bipartite),
          as.integer(Clist$nterms),
          as.character(Clist$fnamestring),
          as.character(Clist$snamestring),
          as.character(MHproposal$name), as.character(MHproposal$package),
          as.double(Clist$inputs), eta=as.double(rep(eta0,MCMCparams$samplesize)),
          as.integer(MCMCparams$samplesize),
          as.integer(MCMCparams$interval),
          s = as.double(t(MCMCparams$stats)),
          as.integer(MCMCparams$burnin),
          as.double(MCMCparams$PILA.steplength), as.double(MCMCparams$PILA.gamma),
          as.integer(verbose), as.integer(MHproposal$bd$attribs),
          as.integer(MHproposal$bd$maxout), as.integer(MHproposal$bd$maxin),
          as.integer(MHproposal$bd$minout), as.integer(MHproposal$bd$minin),
          as.integer(MHproposal$bd$condAllDegExact), as.integer(length(MHproposal$bd$attribs)),
          as.integer(MCMCparams$Clist.miss$heads), as.integer(MCMCparams$Clist.miss$tails),
          as.integer(MCMCparams$Clist.miss$nedges),
          theta.mean.save=if(debug)double((MCMCparams$samplesize+1)*(Clist$nterms+1)),
          XtX.save=if(debug)double((MCMCparams$samplesize+1)*(Clist$nterms+1)^2),
          XtY.save=if(debug)double((MCMCparams$samplesize+1)*(Clist$nterms+1)*Clist$nterms),
          beta.save=if(debug)double((MCMCparams$samplesize+1)*(Clist$nterms+1)*Clist$nterms),
          direction.save=if(debug)double((MCMCparams$samplesize+1)*Clist$nterms),
          dtheta.save=if(debug)double((MCMCparams$samplesize+1)*Clist$nterms),
          insensitive.save=if(debug)integer((MCMCparams$samplesize+1)*Clist$nterms),
          ineffectual.save=if(debug)integer((MCMCparams$samplesize+1)*Clist$nterms),
          dropped.save=if(debug)integer((MCMCparams$samplesize+1)),
          PACKAGE="ergm")
  ## save the results
  copy.named(z)
}
