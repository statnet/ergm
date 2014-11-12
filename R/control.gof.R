#  File R/control.gof.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
#=================================================================
# This file contains the 2 following functions for controlling
# goodness-of-fit computations
#            <control.gof.ergm>
#            <control.gof.formula>
#=================================================================



########################################################################
# Both of the <control.gof.X> functions return a control list for
# customing the fitting procedure used by the gof code
#
# --PARAMETERS--
#   prop.weights  : specifies the method used to allocate probabilities
#                  of being proposed to dyads; options are "TNT",
#                   "random", "nonobserved" and "default"; default=
#                   NULL if X is an ergm (which then uses the weights
#                   that the ergm was fit by); default="default" if
#                   X is a formula (which picks a reasonable default
#                   considering any constraints)
#
# --IGNORED--
#   prop.args     : an alternative, direct way of specifying additional
#                   arguments to proposal; as far as I can tell, the
#                   only use for 'prop.args' is to supply the name
#                   of a nodal attribute for use in the
#                   <InitMHP.nobetweengroupties> function, but this
#                   function will never be called in the path from
#                   <ergm.gof> which is the only code using this
#                   control list.
#   maxchanges    : ??; default=1000000
#
# --RETURNED--
#   a list of the above parameters
#
#########################################################################

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

  set.control.class()
}


control.gof.formula<-function(nsim=100,
                              MCMC.burnin=10000,
                              MCMC.interval=1000,
                              MCMC.prop.weights="default",
                              MCMC.prop.args=list(),
                              
                              MCMC.init.maxedges=20000,
                              MCMC.packagenames=c(),
                              
                              MCMC.runtime.traceplot=FALSE,          
                              network.output="network",
                                                     
                              seed=NULL,
                              parallel=0,
                              parallel.type=NULL,
                              parallel.version.check=TRUE){
  control<-list()
  for(arg in names(formals(sys.function())))
    control[arg]<-list(get(arg))
  
  set.control.class()
}
