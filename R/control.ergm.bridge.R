control.ergm.bridge<-function(nsteps=20, # Number of geometric bridges to use
                              MCMC.burnin=10000,
                              MCMC.interval=100,
                              MCMC.samplesize=10000, # Total number of MCMC draws to use (to be divided up among the bridges, so each bridge gets \code{sample.size/nsteps} draws.
                              obs.MCMC.samplesize=MCMC.samplesize,
                              obs.MCMC.interval=MCMC.interval,
                              obs.MCMC.burnin=MCMC.burnin,
                              
                              MCMC.prop.weights="default",
                              MCMC.prop.args=list(),

                              MCMC.init.maxedges=20000,
                              MCMC.packagenames="ergm",
                              
                              seed=NULL){

  control<-list()
  formal.args<-formals(sys.function())
  for(arg in names(formal.args))
    control[arg]<-list(get(arg))

  set.control.class()
}


