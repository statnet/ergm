stergm.RM <- function(theta.form0, nw, model.form, model.diss, Clist, 
                            theta.diss,
                            MCMCparams, MHproposal.form, MHproposal.diss,
                            verbose=FALSE){
  # This is based on Snijders (2002), J of Social Structure
  # and Snijders and van Duijn (2002) from A Festscrift for Ove Frank
  # Both papers are available from Tom Snijders' web page: 
  #          http://stat.gamma.rug.nl/snijders/publ.htm


  cat("Robbins-Monro algorithm with theta_F_0 = (",theta.form0, ") and theta_D = (",theta.diss,")\n" )
  eta.form0 <- ergm.eta(theta.form0, model.form$etamap)
  eta.diss <- ergm.eta(theta.diss, model.diss$etamap)


  z <- stergm.phase12.C(nw, Clist$meanstats, model.form, model.diss, MHproposal.form, MHproposal.diss,
                        eta.form0, eta.diss, MCMCparams, verbose=verbose)
  ## Phase 3 doesn't give us anything at the moment...
  #MCMCparams$samplesize<-MCMCparams$RM.phase3n
  # Skip burnin this time: if we haven't converged yet, there ain't
  # much anyone can do.
  #MCMCparams$RM.burnin->MCMCparams$time.burnin
  #MCMCparams$RM.interval->MCMCparams$time.interval
  #s <- stergm.getMCMCsample(nw, model.form, model.diss, 
  #                           MHproposal.form, MHproposal.diss, z$eta.form, theta.diss,
  #                           MCMCparams, verbose)
  
  #ve<-with(z,list(coef=eta,sample=s$statsmatrix.form,sample.miss=NULL))
  ve<-with(z,list(coef.form=eta.form,coef.diss=theta.diss))
  
  #endrun <- MCMCparams$burnin+MCMCparams$interval*(ve$samplesize-1)
  #attr(ve$sample, "mcpar") <- c(MCMCparams$burnin+1, endrun, MCMCparams$interval)
  #attr(ve$sample, "class") <- "mcmc"
  
  structure(c(ve, list(newnetwork=nw, 
                       theta.form.original=theta.form0,
                       #interval=MCMCparams$interval, burnin=MCMCparams$burnin, 
                       network=nw)),
            class="stergm")
}

stergm.phase12.C <- function(g, meanstats, model.form, model.diss, 
                             MHproposal.form, MHproposal.diss, eta.form0, eta.diss,
                             MCMCparams, verbose) {
  # ms <- MCMCparams$meanstats
  # if(!is.null(ms)) {
  #   if (is.null(names(ms)) && length(ms) == length(model.form$coef.names))
  #     names(ms) <- model.form$coef.names
  #   obs <- MCMCparams$orig.obs
  #   obs <- obs[match(names(ms), names(obs))]
  #   ms  <-  ms[match(names(obs), names(ms))]
  #   matchcols <- match(names(ms), names(obs))
  #   if (any(!is.na(matchcols))) {
  #     ms[!is.na(matchcols)] <- ms[!is.na(matchcols)] - obs[matchcols[!is.na(matchcols)]]
  #   }
  # }
  Clist.form <- ergm.Cprepare(g, model.form)
  Clist.diss <- ergm.Cprepare(g, model.diss)
  maxchanges <- max(MCMCparams$maxchanges, Clist.form$nedges)/5
  MCMCparams$maxchanges <- MCMCparams$maxchanges/5
  if(verbose){cat(paste("MCMCDyn workspace is",maxchanges,"\n"))}
  
  z <- .C("MCMCDynPhase12",
          # Observed/starting network. 1
          as.integer(Clist.form$heads), as.integer(Clist.form$tails), 
          as.integer(Clist.form$nedges), as.integer(Clist.form$maxpossibleedges),
          as.integer(Clist.form$n),
          as.integer(Clist.form$dir), as.integer(Clist.form$bipartite),
          # Order code. 8
          as.integer(Clist.diss$stergm.order.code),
          # Formation terms and proposals. 9
          as.integer(Clist.form$nterms), as.character(Clist.form$fnamestring), as.character(Clist.form$snamestring),
          as.character(MHproposal.form$name), as.character(MHproposal.form$package),
          as.double(Clist.form$inputs), eta.form=as.double(eta.form0),
          # Formation parameter fitting. 16
          as.double(summary(model.form$formula)-meanstats),
          as.double(MCMCparams$RM.init_gain),
          as.integer(MCMCparams$RM.phase1n_base),
          as.integer(MCMCparams$RM.phase2n_base),
          as.integer(MCMCparams$RM.phase2sub),              
          # Dissolution terms and proposals. 21
          as.integer(Clist.diss$nterms), as.character(Clist.diss$fnamestring), as.character(Clist.diss$snamestring),
          as.character(MHproposal.diss$name), as.character(MHproposal.diss$package),
          as.double(Clist.diss$inputs), as.double(eta.diss),
          # Degree bounds.
          as.integer(MHproposal.form$bd$attribs), 
          as.integer(MHproposal.form$bd$maxout), as.integer(MHproposal.form$bd$maxin),
          as.integer(MHproposal.form$bd$minout), as.integer(MHproposal.form$bd$minin),
          as.integer(MHproposal.form$bd$condAllDegExact), as.integer(length(MHproposal.form$bd$attribs)), 
          # MCMC settings.              
          as.integer(MCMCparams$RM.burnin),
          as.integer(MCMCparams$RM.interval),
          as.integer(MCMCparams$MH.burnin),
          # Space for output.
          as.integer(maxchanges),
          # Verbosity.
          as.integer(verbose), 
          PACKAGE="ergm") 

  eta.form <- z$eta
  names(eta.form) <- names(eta.form0)

  list(meanstats=MCMCparams$meanstats,
       eta.form=eta.form)
}
