#  File R/control.logLik.ergm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
control.logLik.ergm<-function(nsteps=20,
                              MCMC.burnin=NULL,
                              MCMC.interval=NULL,
                              MCMC.samplesize=NULL,
                              obs.MCMC.samplesize=MCMC.samplesize,
                              obs.MCMC.interval=MCMC.interval,
                              obs.MCMC.burnin=MCMC.burnin,
                              
                              MCMC.prop.weights=NULL,
                              MCMC.prop.args=NULL,

                              warn.dyads=TRUE,

                              MCMC.init.maxedges=NULL,
                              MCMC.packagenames=NULL,
                              
                              seed=NULL){

  control<-list()
  formal.args<-formals(sys.function())
  for(arg in names(formal.args))
    control[arg]<-list(get(arg))

  set.control.class()
}
