ergm.robmon.dyn <- function(theta0, nw, model.form, model.diss, Clist, 
                            gamma0,
                            MCMCparams, MHproposal.form, MHproposal.diss,
                            verbose=FALSE){
  # This is based on Snijders (2002), J of Social Structure
  # and Snijders and van Duijn (2002) from A Festscrift for Ove Frank
  # Both papers are available from Tom Snijders' web page: 
  #          http://stat.gamma.rug.nl/snijders/publ.htm

  eta0 <- ergm.eta(theta0, model.form$etamap)

  cat("Robbins-Monro algorithm with theta_0 equal to:\n")
  print(theta0)
  eta0 <- ergm.eta(theta0, model.form$etamap)

  z <- ergm.phase12.dyn(nw, Clist$meanstats, model.form, model.diss, MHproposal.form, MHproposal.diss,
                        eta0, gamma0, MCMCparams, verbose=verbose)

  MCMCparams$samplesize<-MCMCparams$RobMon.phase3n
  # Skip burnin this time: if we haven't converged yet, there ain't
  # much anyone can do.
  MCMCparams$burnin<-0
  s <- ergm.getMCMCDynsample(nw, model.form, model.diss, 
                             MHproposal.form, MHproposal.diss, z$eta, gamma0,
                             MCMCparams, verbose)
  
  ve<-with(z,list(coef=eta,sample=s$statsmatrix.form,sample.miss=NULL))
  
  endrun <- MCMCparams$burnin+MCMCparams$interval*(ve$samplesize-1)
  attr(ve$sample, "mcpar") <- c(MCMCparams$burnin+1, endrun, MCMCparams$interval)
  attr(ve$sample, "class") <- "mcmc"
  
  structure(c(ve, list(newnetwork=nw, 
                       theta.original=theta0,
                       interval=MCMCparams$interval, burnin=MCMCparams$burnin, 
                       network=nw)),
            class="ergm")
}
