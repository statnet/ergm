#  File ergm/R/control.simulate.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
########################################################################
# The <control.simulate.X> functions each create a list of paramaters
# for customizing simulation rountines
########################################################################
control.simulate<-control.simulate.formula<-control.simulate.formula.ergm<-function(MCMC.burnin=1000,
                                                     MCMC.interval=1000,
                                                     MCMC.prop.weights="default",
                                                     MCMC.prop.args=list(),

                                                     MCMC.init.maxedges=20000,
                                                     MCMC.packagenames="ergm",

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

  
  control
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
  
  control
}

control.simulate.ergm.toplevel<-function(control,...){
  ergm.simulate.args<-list(...)
  old.controls<-list(burnin="MCMC.burnin",MCMCsamplesize="MCMLE.samplesize",interval="MCMC.interval")
  for(arg in names(old.controls))
    if(arg %in% names(ergm.simulate.args)){
      warning("Passing ",arg," to simulate.ergm(...) is deprecated and may be removed in a future version. Specify it as control.simulate.ergm(",old.controls[[arg]],"=...) instead.")
      control[old.controls[[arg]]]<-list(ergm.simulate.args[[arg]])
    }
  
  control
}
