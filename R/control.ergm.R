#  File R/control.ergm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
###########################################################################
# The <control.ergm> function allows the ergm fitting process to be tuned
# by returning a list of several control parameters
#
# --PARAMETERS--
#   prop.weights     : the method to allocate probabilities of being proposed
#                      to dyads as "TNT", "random", "nonobserved", or "default"
#                      default="default", which is based upon the ergm constraints
#   prop.args        : a direct way of specifying additional arguments to proposal
#   nr.maxit         : the maximum number of Newton-Raphson iterations to use;
#                      default=100
#   nr.reltol        : the relative tolerance passed to the native R routine
#                      <optim>; default=sqrt(.Machine$double.eps)
#   calc.mcmc.se     : whether MCMC sampling error should be included into the
#                      se estimation (T or F); default=TRUE
#   hessian          : whether the hessian matrix should be computed (T or F);
#                      default=TRUE
#   compress         : whether the stats matrix should be compressed to the set
#                      of unique statistics with a column of probability weights
#                      post-pended; default=FALSE
#   SAN.maxit        : the maximum number of re-runs used to create the 
#                      SAN-ed network and formula; default=10
#   SAN.burnin       : the burnin value used to create the SAN-ed network and
#                      formula; default=NULL
#   maxedges         : the maximum number of edges to allocate space for; default=20000
#   maxchanges       : the maximum number of changes in dynamic network simulation for
#                      which to allocate space; default=1000000
#   obs.MCMCsamplesize:            ; default= NULL,
#   obs.interval    :              ; default= NULL,
#   obs.burnin      :              ; default= NULL,
#   MPLEtype         : the method for MPL estimation as "penalized", "glm" or
#                      "logitreg"; default="glm"
#   steplength       : a multiplier for step length to make fitting more stable at the
#                      cost of efficiency; default=0.5
#   adaptive.trustregion:           ; default=3
#   adaptive.epsilon :              ; default=0.01
#   sequential       : whether the next iteration of the fit should use the last network
#                      sampled as the starting point; the alternative is to always begin
#                      from the orginial network; default=TRUE
#   drop             : whether degenerate terms should be dropped from the fit (T or F);
#                      default=TRUE
#   force.mcmc       : whether ML estimation should be used (T or F); ignored if the model
#                      is not dyadic independent; default=FALSE
#   check.degeneracy : whether the diagnostics should include a degeneacy check (T or F);
#                      default=FALSE
#   mcmc.precision   : a vector of the upper bounds on the precision of the standard errors
#                      induced by the MCMC algorithm; default=0.05
#   metric           : the name of the optimization metric to use, as one of
#                      "Median.Likelihood", "lognormal", "logtaylor","EF.Likelihood" or "naive";
#                      default="Median.Likelihood"
#   method           : the name of the optimaztion method to use, as either "BFGS" or
#                      "Nelder-Mead"; this is an <optim> param; default="BFGS"
#   trustregion      : the maximum amount that the likelihood will be allowed to increase
#                      in an iteration; default=0.5 if 'style'="Stochastic-Approximation",
#                      default=20 otherwise
#   dampening        : (logical) should likelihood dampening be used?
#   dampening.min.ess: effective sample size below which dampening is used
#   dampening.level  : proportional distance from boundary of the convex hull
#                      move
#   style            : the style of ML estimation to use, as one of "Newton-Raphson",
#                      "Robbins-Monro", "Stochastic-Approximation","Stepping";
#                      default="Robbins-Monro"
#   phase1_n         : the number of MCMC samples to draw in Phase 1 of the stochastic
#                      approximation algorithm; default=7 + 3*(# of model terms) if
#                      relevant, otherwise NULL
#   initial_gain     : the initial gain to Phase 2 of the stochastic approximation;
#                      default=0.1 if relevant, otherwise NULL
#   nsubphases       : the number of sub-phases in Phase 2 of the stochastic approximation;
#                      default='maxit'
#   niterations      : the number of MCMC samples to draw in Phase 2 of the stochastic
#                      approximation; default=7 + (# model terms) if relevant, otherwise NULL
#   phase3_n         : the sample size for Phase 3 of the stocastic approximation;
#                      default=1000 if relevant, otherwise NULL
#   RM.init_gain     : this is only used to adjust 'aDdiaginv'in phase1,
#                      in particular:
#                             aDdiaginv = gain/sqrt(aDdiaginv)
#                      default=.5
#   RM.phase1n_base  : this helps define the 'phase1n' param, which in turn
#                      multiplies 'RM.interval' to control the number of
#                      phase1 iterations; this is the base portion of 'phase1n',
#                      which is added to 3*(the number of formation coefficients)
#                      to form 'phase1n'; default=7
#   RM.phase2sub     : phase2 is a 3-deep nested for-loop and 'RM.phase2sub' limits
#                      the outer loop counter; default=7
#   RM.phase2n_base  : this helps define the 'phase2n' param, which in turn
#                      limits the phase2 middle loop counter; this is the
#                      base portion of 'phase2n', which is added to 7+(the number
#                      of formation coefficients) to form 'phase2n'; default=100
#   RM.phase3n   :  ??, i couldn't find this used anywhere; default=500
#   stepMCMCsize       : MCMC sample size for the preliminary steps of the "Stepping"
#                        optimization style; default=100
#   gridsize           : a integer N such that the "Stepping" style of optimization
#                        chooses a step length equal to the largest possible multiple
#                        of 1/N;  default=100
#   packagenames       : the packages in which change statistics are found; default="ergm"
#   parallel           : the number of threads in which to run sampling; default=0
#   returnMCMCstats    : whether the matrix of change stats from the MCMC should be returned as
#                        the mcmc object 'sample'; default=TRUE
#   burnin.retry       : maximum number of times to retry burning in before giving up
#   burnin.check.last  : last what fraction of burnin to check for trending
#   burnin.check.alpha : the alpha for the test
#
# --RETURNED--
#   a list of the above parameters
#
######################################################################################################

control.ergm<-function(drop=TRUE,

                       init=NULL,
                       init.method=NULL,
                       
                       main.method=c("MCMLE","Robbins-Monro",
                               "Stochastic-Approximation","Stepping"),
                       force.main=FALSE,
                       main.hessian=TRUE,

                       MPLE.max.dyad.types=1e+6, 
                       MPLE.samplesize=50000,                       
                       MPLE.type=c("glm", "penalized"),
                      
                       MCMC.prop.weights="default", MCMC.prop.args=list(),
                       MCMC.interval=1024,
                       MCMC.burnin=MCMC.interval*16,
                       MCMC.samplesize=1024,
                       MCMC.effectiveSize=NULL,
                       MCMC.effectiveSize.damp=10,
                       MCMC.effectiveSize.maxruns=1000,
                       MCMC.effectiveSize.base=1/2,
                       MCMC.effectiveSize.points=5,
                       MCMC.effectiveSize.order=1,
                       MCMC.return.stats=TRUE,
                       MCMC.runtime.traceplot=FALSE,
                       MCMC.init.maxedges=20000,
                       MCMC.max.maxedges=Inf,
                       MCMC.addto.se=TRUE,
                       MCMC.compress=FALSE,
                       MCMC.packagenames=c(),

                       SAN.maxit=10,
                       SAN.burnin.times=10,
                       SAN.control=control.san(coef=init,
                         SAN.prop.weights=MCMC.prop.weights,
                         SAN.prop.args=MCMC.prop.args,
                         SAN.init.maxedges=MCMC.init.maxedges,
                         
                         SAN.burnin=MCMC.burnin*SAN.burnin.times,
                         SAN.interval=MCMC.interval,
                         SAN.packagenames=MCMC.packagenames,
                         MPLE.max.dyad.types=MPLE.max.dyad.types,

                         parallel=parallel,
                         parallel.type=parallel.type,
                         parallel.version.check=parallel.version.check),
                       
                       MCMLE.termination=c("Hummel", "Hotelling", "precision", "none"),
                       MCMLE.maxit=20,
                       MCMLE.conv.min.pval=0.5,
                       MCMLE.NR.maxit=100,
                       MCMLE.NR.reltol=sqrt(.Machine$double.eps),
                       obs.MCMC.samplesize=MCMC.samplesize,
                       obs.MCMC.interval=MCMC.interval,
                       obs.MCMC.burnin=MCMC.burnin,
                       obs.MCMC.burnin.min=obs.MCMC.burnin/10,
                       obs.MCMC.prop.weights=MCMC.prop.weights, obs.MCMC.prop.args=MCMC.prop.args,

                       MCMLE.check.degeneracy=FALSE,
                       MCMLE.MCMC.precision=0.005,
                       MCMLE.MCMC.max.ESS.frac=0.1,
                       MCMLE.metric=c("lognormal", "logtaylor",
                         "Median.Likelihood",
                         "EF.Likelihood", "naive"),
                       MCMLE.method=c("BFGS","Nelder-Mead"),
                       MCMLE.trustregion=20,
                       MCMLE.dampening=FALSE,
                       MCMLE.dampening.min.ess=20,
                       MCMLE.dampening.level=0.1,
                       MCMLE.steplength.margin=0.05,
                       MCMLE.steplength=if(is.null(MCMLE.steplength.margin)) 0.5 else 1,
                       MCMLE.adaptive.trustregion=3,
                       MCMLE.sequential=TRUE,
                       MCMLE.density.guard.min=10000,
                       MCMLE.density.guard=exp(3),
                       MCMLE.effectiveSize=NULL,
                       MCMLE.last.boost=4,
                       MCMLE.Hummel.esteq=TRUE, 
                       MCMLE.Hummel.miss.sample=100,
                       MCMLE.Hummel.maxit=25, 
                       MCMLE.steplength.min=0.0001,
                       
                       SA.phase1_n=NULL, SA.initial_gain=NULL, 
                       SA.nsubphases=4,
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

                       CD.nsteps=8,
                       CD.multiplicity=1,
                       CD.nsteps.obs=128,
                       CD.multiplicity.obs=1,
                       CD.maxit=60,
                       CD.conv.min.pval=0.5,
                       CD.NR.maxit=100,
                       CD.NR.reltol=sqrt(.Machine$double.eps),
                       CD.metric=c("naive", "lognormal", "logtaylor",
                         "Median.Likelihood",
                         "EF.Likelihood"),
                       CD.method=c("BFGS","Nelder-Mead"),
                       CD.trustregion=20,
                       CD.dampening=FALSE,
                       CD.dampening.min.ess=20,
                       CD.dampening.level=0.1,
                       CD.steplength.margin=0.5,
                       CD.steplength=1,
                       CD.adaptive.trustregion=3,
                       CD.adaptive.epsilon=0.01,
                       CD.Hummel.esteq=TRUE, 
                       CD.Hummel.miss.sample=100,
                       CD.Hummel.maxit=25, 
                       CD.steplength.min=0.0001,
                       
                       loglik.control=control.logLik.ergm(),

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
                       obs.interval="obs.MCMC.interval",
                       obs.burnin="obs.MCMC.burnin",
                       compress="MCMC.compress",
                       metric="MCMLE.metric",
                       force.mcmc="force.main",
                       adaptive.trustregion="MCMLE.adaptive.trustregion",
                       adaptive.epsilon="MCMLE.adaptive.epsilon",
                       mcmc.precision="MCMLE.MCMC.precision",
                       method="MCMLE.method",
                       MPLEtype="MPLE.type",
                       check.degeneracy="MCMLE.check.degeneracy",
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

  match.arg.pars <- c("MPLE.type","MCMLE.metric","MCMLE.method","main.method",'MCMLE.termination',"CD.metric","CD.method")
  
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

  if((MCMLE.steplength!=1 || is.null(MCMLE.steplength.margin)) && MCMLE.termination %in% c("Hummel", "precision"))
    stop("Hummel and precision-based termination require non-null MCMLE.steplength.margin and MCMLE.steplength = 1.")
  
  set.control.class()
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
