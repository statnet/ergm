#  File ergm/R/control.logLik.ergm.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
control.logLik.ergm<-function(nsteps=20,
                              MCMC.burnin=NULL,
                              MCMC.interval=NULL,
                              MCMC.samplesize=NULL,
                              obs.MCMC.samplesize=MCMC.samplesize,
                              obs.MCMC.interval=MCMC.interval,
                              obs.MCMC.burnin=MCMC.burnin,
                              
                              MCMC.prop.weights=NULL,
                              MCMC.prop.args=NULL,

                              MCMC.init.maxedges=NULL,
                              MCMC.packagenames=NULL,
                              
                              seed=NULL){

  control<-list()
  formal.args<-formals(sys.function())
  for(arg in names(formal.args))
    control[[arg]]<-get(arg)

  control
}
