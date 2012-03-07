#  File ergm/R/control.gof.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
########################################################################
# Both of the <control.gof.X> functions return a control list for
# customing the fitting procedure used by the gof code
########################################################################
control.gof.ergm<-function(nsim=100,
                           MCMC.burnin=NULL,
                           MCMC.interval=NULL,
                           MCMC.prop.weights=NULL,
                           MCMC.prop.args=NULL,
                           
                           MCMC.init.maxedges=NULL,
                           MCMC.packagenames=NULL,

                           MCMC.runtime.traceplot=FALSE,
                           network.output="network",

                           seed=NULL,
                           parallel=0,
                           parallel.type=NULL,
                           parallel.version.check=TRUE){
  control<-list()
  for(arg in names(formals(sys.function())))
    control[arg]<-list(get(arg))
  control
}


control.gof.formula<-function(nsim=100,
                              MCMC.burnin=1000,
                              MCMC.interval=1000,
                              MCMC.prop.weights="default",
                              MCMC.prop.args=list(),
                              
                              MCMC.init.maxedges=20000,
                              MCMC.packagenames="ergm",
                              
                              MCMC.runtime.traceplot=FALSE,          
                              network.output="network",
                                                     
                              seed=NULL,
                              parallel=0,
                              parallel.type=NULL,
                              parallel.version.check=TRUE){
  control<-list()
  for(arg in names(formals(sys.function())))
    control[arg]<-list(get(arg))
  control
}
