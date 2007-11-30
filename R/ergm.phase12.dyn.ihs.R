ergm.phase12.dyn <- function(g, meanstats, model.form, model.diss, 
                             MHproposal.form, MHproposal.diss, eta0, gamma0,
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
          # Observed/starting network.
          as.integer(Clist.form$heads), as.integer(Clist.form$tails), 
          as.integer(Clist.form$nedges), as.integer(Clist.form$n),
          as.integer(Clist.form$dir), as.integer(Clist.form$bipartite),
          # Order code.
          as.integer(Clist.diss$order.code),
          # Formation terms and proposals.
          as.integer(Clist.form$nterms), as.character(Clist.form$fnamestring), as.character(Clist.form$snamestring),
          as.character(MHproposal.form$name), as.character(MHproposal.form$package),
          as.double(Clist.form$inputs), eta=as.double(eta0),
          # Formation parameter fitting.
          as.double(summary(model.form$formula)-meanstats),
          as.double(MCMCparams$RobMon.init_gain),
          as.integer(MCMCparams$RobMon.phase1n_base),
          as.integer(MCMCparams$RobMon.phase2n_base),
          as.integer(MCMCparams$RobMon.phase2sub),              
          # Dissolution terms and proposals.
          as.integer(Clist.diss$nterms), as.character(Clist.diss$fnamestring), as.character(Clist.diss$snamestring),
          as.character(MHproposal.diss$name), as.character(MHproposal.diss$package),
          as.double(Clist.diss$inputs), as.double(gamma0),
          # Degree bounds.
          as.integer(MHproposal.form$bd$attribs), 
          as.integer(MHproposal.form$bd$maxout), as.integer(MHproposal.form$bd$maxin),
          as.integer(MHproposal.form$bd$minout), as.integer(MHproposal.form$bd$minin),
          as.integer(MHproposal.form$bd$condAllDegExact), as.integer(length(MHproposal.form$bd$attribs)), 
          # MCMC settings.              
          as.integer(MCMCparams$burnin), as.integer(MCMCparams$dyninterval),
          as.integer(MCMCparams$interval),
          # Space for output.
          as.integer(maxchanges),
          # Verbosity.
          as.integer(verbose), 
          PACKAGE="ergm") 

  eta <- z$eta
  names(eta) <- names(eta0)

  list(meanstats=MCMCparams$meanstats,
       eta=eta)
}
