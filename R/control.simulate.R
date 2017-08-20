#  File R/control.simulate.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
#==========================================================
# This file contains the following 2 functions for
# controlling simulation routines
#        <control.simulate> = <control.simulate.formula>
#        <control.simulate.ergm>
#==========================================================


########################################################################
# The <control.simulate.X> functions each create a list of paramaters
# for customizing simulation rountines
#
# --PARAMETERS--
#   prop.weights  : specifies the method used to allocate probabilities
#                   of being proposed to dyads; options are "TNT",
#                   "random", "nonobserved" and "default"; default=
#                   NULL if X is an ergm (which then uses the weights
#                   that the ergm was fit by); default="default" if
#                   X is a formula (which picks a reasonable default
#                   considering any constraints)
#   maxedges      : the maximum number of edges expected in network;
#                   default=20000
#   packagenames  : the names of packages to load on the created
#                   cluster when using parallel threads; currently,
#                   the only recognized package name is "ergm";
#                   default="ergm"
#   network.output: the R class with which to output networks; the
#                   options are "NULL", "network" and "edgelist.compressed"
#                   (which saves space but only supports networks
#                   without vertex attributes); "NULL" does not
#                   return networks; default="network"
#   parallel      : number of threads in which to run the sampling
#
# --IGNORED--
#   prop.args     : an alternative, direct way of specifying additional
#                   arguments to proposal; as far as I can tell, the
#                   only use for 'prop.args' is to supply the name
#                   of a nodal attribute for use in the
#                   <InitMHP.nobetweengroupties> function, but this
#                   function will never be called in the path from
#                   <ergm.san> which is the only code using this
#                   control list.
#   maxchanges    : ??; default=1000000 
#
#
# --RETURNED--
#   a list of the above parameters
#
#########################################################################

control.simulate<-control.simulate.formula<-control.simulate.formula.ergm<-function(MCMC.burnin=10000,
                                                     MCMC.interval=1000,
                                                     MCMC.prop.weights="default",
                                                     MCMC.prop.args=list(),

                                                     MCMC.init.maxedges=20000,
                                                     MCMC.packagenames=c(),

                                                     MCMC.runtime.traceplot=FALSE,  
                                                     network.output="network",
                                                     
                                                     parallel=0,
                                                     parallel.type=NULL,
                                                     parallel.version.check=TRUE,
                                                     ...){
  old.controls <- list(
                       maxedges="MCMC.init.maxedges",
                       prop.weights="MCMC.prop.weights",
                       prop.args="MCMC.prop.args",
                       packagenames="MCMC.packagenames"
                       )

  control<-list()
  formal.args<-formals(sys.function())
  formal.args[["..."]]<-NULL
  for(arg in names(formal.args))
    control[arg]<-list(get(arg))

  for(arg in names(list(...))){
    if(!is.null(old.controls[[arg]])){
      warning("Passing ",arg," to control.simulate.formula(...) is deprecated and may be removed in a future version. Specify it as control.simulate.formula(",old.controls[[arg]],"=...) instead.")
      control[old.controls[[arg]]]<-list(list(...)[[arg]])
    }else{
      stop("Unrecognized control parameter: ",arg,".")
    }
  }

  set.control.class("control.simulate.formula")
}

control.simulate.ergm<-function(MCMC.burnin=NULL,
                                MCMC.interval=NULL,
                                MCMC.prop.weights=NULL,
                                MCMC.prop.args=NULL,

                                MCMC.init.maxedges=NULL,
                                MCMC.packagenames=NULL,
                                
                                MCMC.runtime.traceplot=FALSE,
                                network.output="network",

                                parallel=0,
                                parallel.type=NULL,
                                parallel.version.check=TRUE,
                                ...){
  old.controls <- list(
                       maxedges="MCMC.init.maxedges",
                       prop.weights="MCMC.prop.weights",
                       prop.args="MCMC.prop.args",
                       packagenames="MCMC.packagenames"
                       )

  control<-list()
  formal.args<-formals(sys.function())
  formal.args[["..."]]<-NULL
  for(arg in names(formal.args))
    control[arg]<-list(get(arg))

  for(arg in names(list(...)))
    if(!is.null(old.controls[[arg]])){
      warning("Passing ",arg," to control.simulate.ergm(...) is deprecated and may be removed in a future version. Specify it as control.simulate.ergm(",old.controls[[arg]],"=...) instead.")
      control[old.controls[[arg]]]<-list(list(...)[[arg]])
    }
 
  set.control.class("control.simulate.ergm")
}
