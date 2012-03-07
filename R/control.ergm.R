#  File ergm/R/control.ergm.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
###########################################################################
# The <control.ergm> function allows the ergm fitting process to be tuned
# by returning a list of several control parameters
###########################################################################
control.ergm<-function(drop=TRUE,

                       init=NULL,
                       init.method=NULL,
                       
                       main.method=c("MCMLE","Robbins-Monro",
                               "Stochastic-Approximation","Stepping"),
                       force.main=FALSE,
                       main.hessian=TRUE,

                       MPLE.max.dyad.types=1e+6, 
                       MPLE.max.samplesize=100000,
                       MPLE.samplesize=50000,                       
                       MPLE.type=c("glm", "penalized"),
                      
                       MCMC.prop.weights="default", MCMC.prop.args=list(),
                       MCMC.burnin=10000,
                       MCMC.interval=100,
                       MCMC.samplesize=10000,
                       MCMC.return.stats=TRUE,
                       MCMC.burnin.retries=0,
                       MCMC.burnin.check.last=1/2,
                       MCMC.burnin.check.alpha=0.01,
                       MCMC.runtime.traceplot=FALSE,
                       MCMC.init.maxedges=20000,
                       MCMC.addto.se=TRUE,
                       MCMC.compress=FALSE,
                       MCMC.packagenames="ergm",

                       SAN.maxit=10,
                       SAN.control=control.san(coef=init,
                         SAN.prop.weights=MCMC.prop.weights,
                         SAN.prop.args=MCMC.prop.args,
                         SAN.init.maxedges=MCMC.init.maxedges,
                         
                         SAN.burnin=MCMC.burnin,
                         SAN.interval=MCMC.interval,
                         SAN.packagenames=MCMC.packagenames,

                         parallel=parallel,
                         parallel.type=parallel.type,
                         parallel.version.check=parallel.version.check),

                       MCMLE.maxit=20,
                       MCMLE.conv.min.pval=0.5,
                       MCMLE.NR.maxit=100,
                       MCMLE.NR.reltol=sqrt(.Machine$double.eps),
                       MCMLE.obs.MCMC.samplesize=MCMC.samplesize,
                       MCMLE.obs.MCMC.interval=MCMC.interval,
                       MCMLE.obs.MCMC.burnin=MCMC.burnin,
                       MCMLE.check.degeneracy=FALSE,
                       MCMLE.MCMC.precision=0.05,
                       MCMLE.metric=c("lognormal", "Median.Likelihood",
                         "EF.Likelihood", "naive"),
                       MCMLE.method=c("BFGS","Nelder-Mead"),
                       MCMLE.trustregion=20,
                       MCMLE.steplength=0.5,
                       MCMLE.adaptive.trustregion=3,
                       MCMLE.adaptive.epsilon=0.01,
                       MCMLE.sequential=TRUE,

                       SA.phase1_n=NULL, SA.initial_gain=NULL, 
                       SA.nsubphases=MCMLE.maxit,
                       SA.niterations=NULL, 
                       SA.phase3_n=NULL,
                       SA.trustregion=0.5,

                       RM.phase1n_base=7,
                       RM.phase2n_base=100,
                       RM.phase2sub=7,
                       RM.init_gain=0.5,
                       RM.phase3n=500,

                       Step.MCMC.samplesize=100,
                       Step.maxit=50,
                       Step.gridsize=100,

                       loglik.control=control.ergm.bridge(MCMC.burnin=MCMC.burnin,
                         MCMC.interval=MCMC.interval,
                         MCMC.samplesize=MCMC.samplesize,
                         MCMC.prop.weights=MCMC.prop.weights,
                         MCMC.prop.args=MCMC.prop.args,
                         MCMC.init.maxedges=MCMC.init.maxedges,
                         MCMC.packagenames=MCMC.packagenames),

                       seed=NULL,
                       parallel=0,
                       parallel.type=NULL,
                       parallel.version.check=TRUE,
                       
                       ...
                       ){
  old.controls <- list(nr.maxit="MCMLE.NR.maxit",
                       nr.reltol="MCMLE.NR.reltol",
                       maxNumDyadTypes="MPLE.max.dyad.types",
                       maxedges="MCMC.init.maxedges",
                       steplength="MCMLE.steplength",
                       initialfit="init.method",
                       style="main.method",
                       obs.MCMCsamplesize="MCMLE.obs.samplesize",
                       obs.interval="MCMLE.obs.MCMC.interval",
                       obs.burnin="MCMLE.obs.MCMC.burnin",
                       compress="MCMC.compress",
                       metric="MCMLE.metric",
                       force.mcmc="force.main",
                       adaptive.trustregion="MCMLE.adaptive.trustregion",
                       adaptive.epsilon="MCMLE.adaptive.epsilon",
                       mcmc.precision="MCMLE.mcmc.precision",
                       method="MCMLE.method",
                       MPLEtype="MPLE.type",
                       check.degeneracy="MCMLE.check.degeneracy",
                       maxMPLEsamplesize="MPLE.max.samplesize",
                       MPLEsamplesize="MPLE.samplesize",
                       phase1_n="SA.phase1_n", initial_gain="SA.initial_gain", 
                       nsubphases="SA.nsubphases", niterations="SA.niterations", phase3_n="SA.phase3_n",
                       RobMon.phase1n_base="RM.phase1n_base",
                       RobMon.phase2n_base="RM.phase2n_base",
                       RobMon.phase2sub="RM.phase2sub",
                       RobMon.init_gain="RM.init_gain",
                       RobMon.phase3n="RM.phase3n",
                       trustregion="MCMLE.trustregion",
                       stepMCMCsize="Step.MCMC.samplesize",
                       steppingmaxit="Step.maxit",
                       gridsize="Step.gridsize",
                       sequential="MCMLE.sequential",
                       returnMCMCstats="MCMC.return.stats",
                       calc.mcmc.se="MCMC.addto.se",
                       hessian="main.hessian",
                       prop.weights="MCMC.prop.weights",
                       prop.args="MCMC.prop.args",
                       packagenames="MCMC.packagenames"
                       )

  match.arg.pars=c("MPLE.type","MCMLE.metric","MCMLE.method","main.method")
  
  control<-list()
  formal.args<-formals(sys.function())
  formal.args[["..."]]<-NULL
  for(arg in names(formal.args))
    control[arg]<-list(get(arg))

  for(arg in names(list(...))){
    if(!is.null(old.controls[[arg]])){
      warning("Passing ",arg," to control.ergm(...) is deprecated and may be removed in a future version. Specify it as control.ergm(",old.controls[[arg]],"=...) instead.")
      control[old.controls[[arg]]]<-list(list(...)[[arg]])
    }else{
      stop("Unrecognized control parameter: ",arg,".")
    }
  }

  for(arg in match.arg.pars)
    control[arg]<-list(match.arg(control[[arg]][1],eval(formal.args[[arg]])))
  
  control
}

control.ergm.toplevel<-function(control,...){
  ergm.args<-list(...)
  old.controls<-list(burnin="MCMC.burnin",MCMCsamplesize="MCMC.samplesize",interval="MCMC.interval",maxit="MCMLE.maxit",seed="seed",theta0="init",coef="init")
  for(arg in names(old.controls))
    if(arg %in% names(ergm.args)){
      warning("Passing ",arg," to ergm(...) is deprecated and may be removed in a future version. Specify it as control.ergm(",old.controls[[arg]],"=...) instead.")
      control[old.controls[[arg]]]<-list(ergm.args[[arg]])
    }
  
  control
}
