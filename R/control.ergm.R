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
#   SAN.burnin       : the burnin value used to create the SAN-ed network and
#                      formula; default=NULL
#   maxNumDyadTypes  : the maximum number of unique psuedolikelihood change stats
#                      to be allowed if 'compress'=TRUE; ignored if 'compress'!=TRUE;
#                      default=1e+6
#   maxedges         : the maximum number of edges to allocate space for; default=20000
#   maxchanges       : the maximum number of changes in dynamic network simulation for
#                      which to allocate space; default=1000000
#   maxMPLEsamplesize: the sample size to use for endogenous sampling in the psuedo-
#                      likelihood computation; default=100000
#   obs.MCMCsamplesize:            ; default= NULL,
#   obs.interval    :              ; default= NULL,
#   obs.burnin      :              ; default= NULL,
#   MPLEtype         : the method for MPL estimation as "penalized", "glm" or
#                      "logitreg"; default="glm"
#   trace            : the number of levels of tracing information to produce during
#                      optimization; see <?optim> for details; default=0
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
#                      "Median.Likelihood", "lognormal", "EF.Likelihood" or "naive";
#                      default="Median.Likelihood"
#   method           : the name of the optimaztion method to use, as either "BFGS" or
#                      "Nelder-Mead"; this is an <optim> param; default="BFGS"
#   trustregion      : the maximum amount that the likelihood will be allowed to increase
#                      in an iteration; default=0.5 if 'style'="Stochastic-Approximation",
#                      default=20 otherwise
#   initial.loglik   : an initial value for the log-likelihood; default=NULL
#   initial.network  : an initial network for the MCMC procedure; default=NULL
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
#   RobMon.phase3n   :  ??, i couldn't find this used anywhere; default=500
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

control.ergm<-function(prop.weights="default",prop.args=list(),
                       nr.maxit=100,
                       nr.reltol=sqrt(.Machine$double.eps),
                       calc.mcmc.se=TRUE, hessian=TRUE,
                       compress=FALSE,
                       SAN.burnin=NULL,
                       maxNumDyadTypes=1e+6, 
                       maxedges=20000,
                       maxchanges=1000000,
                       initialfit=NULL,
                       maxMPLEsamplesize=100000,
                       MPLEsamplesize=50000,
                       obs.MCMCsamplesize=NULL,
                       obs.interval=NULL,
                       obs.burnin=NULL,
                       MPLEtype=c("glm", "penalized"),
                       trace=0,
                       steplength=0.5,
                       adaptive.trustregion=3,
                       adaptive.epsilon=0.01,
                       sequential=TRUE,
                       drop=TRUE,
                       force.mcmc=FALSE,
                       check.degeneracy=FALSE,
                       mcmc.precision=0.05,
                       metric=c("lognormal", "Median.Likelihood",
                                "EF.Likelihood", "naive"),
                       method=c("BFGS","Nelder-Mead"),
                       trustregion=20,
                       loglik.nsteps=20,
                       initial.network=NULL,
                       style=c("Newton-Raphson","Robbins-Monro",
                               "Stochastic-Approximation","Stepping"),
                       phase1_n=NULL, initial_gain=NULL, 
                       nsubphases="maxit", niterations=NULL, phase3_n=NULL,
                       RobMon.phase1n_base=7,
                       RobMon.phase2n_base=100,
                       RobMon.phase2sub=7,
                       RobMon.init_gain=0.5,
                       RobMon.phase3n=500,
                       stepMCMCsize=100,
                       gridsize=100,
                       packagenames="ergm",
                       parallel=0,
                       parallel.type=NULL,
                       returnMCMCstats=TRUE,
                       burnin.retries=0,
                       burnin.check.last=1/2,
                       burnin.check.alpha=0.01,
                       runtime.traceplot=FALSE
                       ){
  control<-list()
  for(arg in names(formals(sys.function())))
    control[[arg]]<-get(arg)
  
  control$MPLEtype<-match.arg(MPLEtype)
  control$metric<-match.arg(metric)
  control$method<-match.arg(method)
  control$style<-match.arg(style)
  control$nsubphases<-match.arg(nsubphases)
  if(missing(trustregion) & control$style=="Stochastic-Approximation"){
   control$trustregion <- 0.5
  }

  control
}
