ergm.runPILAsampler <- function(nw, model, MHproposal, eta0, MCMCparams, 
                               verbose) {

  Clist <- ergm.Cprepare(nw, model)
  z <- ergm.PILAslave(Clist,MHproposal,eta0,MCMCparams,verbose)
  nedges <- z$newnwheads[1]

  statsmatrix <- matrix(z$s, nrow=MCMCparams$samplesize,
                        ncol=Clist$nparam,
                        byrow = TRUE)
  #statsmatrix.b <- matrix(z$s.b, nrow=MCMCparams$samplesize,
  #                      ncol=Clist$nparam,
  #                      byrow = TRUE)
  etamatrix <- matrix(z$eta, nrow=MCMCparams$samplesize,
                      ncol=Clist$nparam,
                      byrow=TRUE)
  #etamatrix.b <- matrix(z$eta.b, nrow=MCMCparams$samplesize,
  #                    ncol=Clist$nparam,
  #                    byrow=TRUE)
  colnames(statsmatrix) <- model$coef.names

  list(statsmatrix=statsmatrix, etamatrix=etamatrix,
       #statsmatrix.b=statsmatrix.b, etamatrix.b=etamatrix.b,
       #newnetwork=newnetwork, 
       meanstats=Clist$meanstats, nedges=nedges)
}
# Function the slaves will call to perform a validation on the
# mcmc equal to their slave number.
# Assumes: Clist MHproposal eta0 MCMCparams maxedges verbose
ergm.PILAslave <- function(Clist,MHproposal,eta0,MCMCparams,verbose) {
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
          double(0), #eta.b=double(MCMCparams$burnin*length(eta0)),
          double(0), #s.b=double(MCMCparams$burnin*length(eta0)),
          as.double(MCMCparams$PILA.steplength), as.double(MCMCparams$PILA.gamma),
          as.integer(verbose), as.integer(MHproposal$bd$attribs),
          as.integer(MHproposal$bd$maxout), as.integer(MHproposal$bd$maxin),
          as.integer(MHproposal$bd$minout), as.integer(MHproposal$bd$minin),
          as.integer(MHproposal$bd$condAllDegExact), as.integer(length(MHproposal$bd$attribs)),
          as.integer(MCMCparams$Clist.miss$heads), as.integer(MCMCparams$Clist.miss$tails),
          as.integer(MCMCparams$Clist.miss$nedges),
          PACKAGE="ergm")
                                        # save the results
  list(s=z$s, eta=z$eta, s.b=z$s.b, eta.b=z$eta.b, newnwheads=z$newnwheads, newnwtails=z$newnwtails)
}
