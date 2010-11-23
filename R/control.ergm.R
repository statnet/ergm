###########################################################################
# The <control.ergm> function allows the ergm fitting process to be tuned
# by returning a list of several control parameters
#
# --PARAMETERS--
#   prop.weights     : the method to allocate probabilities of being proposed
#                      to dyads as "TNT", "random", "nonobserved", or "default"
#                      default="default", which is based upon the ergm constraints
#   prop.args        : an alternative, direct way of specifying additional
#                      arguments to proposal
#   prop.weights.diss: as 'prop.weights', but for the dissolution model
#   prop.args.diss   : as 'prop.args', but for the dissoultion model
#   nr.maxit         : the maximum number of Newton-Raphson iterations to use;
#                      default=100
#   calc.mcmc.se     : whether MCMC sampling error should be included into the
#                      se estimation (T or F); default=TRUE
#   hessian          : whether the hessian matrix should be computed (T or F);
#                      default=TRUE
#   compress         : whether the stats matrix should be compressed to the set
#                      of unique statistics with a column of probability weights 
#                      post-pended; default=TRUE
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
#   MPLEtype         : the method for MPL estimation as "penalized", "glm" or
#                      "logitreg"; default="glm"
#   nr.reltol        : the relative tolerance passed to the native R routine
#                      <optim>; default=sqrt(.Machine$double.eps)
#   trace            : the number of levels of tracing information to produce during
#                      optimization; see <?optim> for details; default=0
#   steplength       : a multiplier for step length to make fitting more stable at the
#                      cost of efficiency; default=0.5
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
#                      in an iteration; default=20
#   initial.loglik   : an initial value for the log-likelihood; default=NULL
#   initial.network  : an initial network for the MCMC procedure; default=NULL
#   style            : the style of ML estimation to use, as one of "Newton-Raphson",
#                      "Robbins-Monro", "Stochastic-Approximation","Stepping" or "PILA";
#                      default="Robbins-Monro"
#   style.dyn        : the style of method of moments estimation to use, as "Robbins-Monro",
#                      "SPSA" or "SPSA2"; "Robbins-Monro" should only be used if it is
#                      known a priori that the derivative of each element of the equilibrium
#                      expected values of statistics of interest with respect to the
#                      corresponding formation phase parameter is positive;
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
#   RobMon.phase1n_base: ??; default=7
#   RobMon.phase2n_base: ??; default=100
#   RobMon.phase2sub   : ??; default=7
#   RobMon.init_gain   : ??; default=0.5
#   RobMon.phase3n     : ??; default=500
#   stepMCMCsize       : MCMC sample size for the preliminary steps of the "Stepping"
#                        optimization style; default=100
#   gridsize           : a integer N such that the "Stepping" style of optimization
#                        chooses a step length equal to the largest possible multiple
#                        of 1/N;  default=100
#   dyninterval        : number of MH proposals for each phase in the dynamic network
#                        simulation; default=1000
#   packagenames       : the packages in which change statistics are found; default="ergm"
#   parallel           : the number of threads in which to run sampling; default=0
#   returnMCMCstats    : whether the matrix of change stats from the MCMC should be returned as
#                        the mcmc object 'sample'; default=TRUE
#
# --RETURNED--
#   a list of the above parameters
#
######################################################################################################

control.ergm<-function(prop.weights="default",prop.args=NULL,
                       prop.weights.diss="default",prop.args.diss=NULL,
                       nr.maxit=100,
                       calc.mcmc.se=TRUE, hessian=TRUE,
                       compress=TRUE,
                       SAN.burnin=NULL,
                       maxNumDyadTypes=1e+6, 
                       maxedges=20000,
                       maxchanges=1000000,
                       maxMPLEsamplesize=100000,
                       MPLEtype=c("glm", "penalized"),
                       nr.reltol=sqrt(.Machine$double.eps),
                       trace=0,
                       steplength=0.5,
                       sequential=TRUE,
                       drop=TRUE,
                       force.mcmc=FALSE,
                       check.degeneracy=FALSE,
                       mcmc.precision=0.05,
                       metric=c("lognormal", "Median.Likelihood",
                                "EF.Likelihood", "naive"),
                       method=c("BFGS","Nelder-Mead"),
                       trustregion=20,
                       initial.loglik=NULL,
                       initial.network=NULL,
                       style=c("Newton-Raphson","Robbins-Monro",
                               "Stochastic-Approximation","Stepping","PILA"),
                       style.dyn=c("Robbins-Monro","SPSA", "SPSA2"),
                       phase1_n=NULL, initial_gain=NULL, 
                       nsubphases="maxit", niterations=NULL, phase3_n=NULL,
                       RobMon.phase1n_base=7,
                       RobMon.phase2n_base=100,
                       RobMon.phase2sub=7,
                       RobMon.init_gain=0.5,
                       RobMon.phase3n=500,
                       stepMCMCsize=100,
                       gridsize=100,
                       dyninterval=1000,
                       packagenames="ergm",
                       parallel=0,
                       returnMCMCstats=TRUE){
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
