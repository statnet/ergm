###########################################################################
# The <control.stergm> function allows the ergm fitting process to be tuned
# by returning a list of several control parameters
#
# --PARAMETERS--
#   prop.weights.form: the method to allocate probabilities of being proposed
#                      to dyads in the formation stage, as "TNT", "random",
#                      "nonobserved", or "default"
#                      default="default", which is based upon the ergm constraints
#   prop.args.form   : an alternative, direct way of specifying additional
#                      arguments to the formation proposals
#   prop.weights.diss: as 'prop.weights.form', but for the dissolution model
#   prop.args.diss   : as 'prop.args.form', but for the dissoultion model
#   compress         : whether the stats matrix should be compressed to the set
#                      of unique statistics with a column of probability weights
#                      post-pended; default=FALSE
#   SAN.burnin       : the burnin value used to create the SAN-ed network and
#                      formula; default=NULL
#   SAN.interval     :
#   maxNumDyadTypes  : the maximum number of unique psuedolikelihood change stats
#                      to be allowed if 'compress'=TRUE; ignored if 'compress'!=TRUE;
#                      default=1e+6
#   maxedges         : the maximum number of edges to allocate space for; default=20000
#   maxchanges       : the maximum number of changes in dynamic network simulation for
#                      which to allocate space; default=1000000
#   maxMPLEsamplesize: the sample size to use for endogenous sampling in the psuedo-
#                      likelihood computation; default=100000
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
#   trace            : the number of levels of tracing information to produce during
#                      optimization; see <?optim> for details; default=0
#   sequential       : whether the next iteration of the fit should use the last network
#                      sampled as the starting point; the alternative is to always begin
#                      from the orginial network; default=TRUE
#   style            : the style of ML estimation to use, as one of "Newton-Raphson",
#                      "Robbins-Monro", "Stochastic-Approximation", or "Stepping";
#                      default="Robbins-Monro"
#   RM.init_gain     : this is only used to adjust 'aDdiaginv'in phase1,
#                      in particular:
#                             aDdiaginv = gain/sqrt(aDdiaginv)
#                      default=0.5
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
#   RM.burnin        : the number of MCMC steps to disregard for the burn-in
#                      period; default=1000
#   RM.interval      : like the SPSA.interval, this seems a little more like
#                      a sample size, than an interval, it helps control the
#                      number of MCMCsteps used in phase1 and phase2; in
#                      phase2, this limits the innermost loop counter; default=100
#   packagenames     : the packages in which change statistics are found; default="ergm"
#   parallel         : the number of threads in which to run sampling; default=0
#
# --RETURNED--
#   a list of the above parameters
#
######################################################################################################

control.stergm<-function(init.form=NULL,
                         init.diss=NULL,
                         init.method=NULL,

                         MCMC.prop.weights.form="default",MCMC.prop.args.form=NULL,
                         MCMC.prop.weights.diss="default",MCMC.prop.args.diss=NULL,
                         MCMC.init.maxedges=20000,
                         MCMC.init.maxchanges=20000,
                         MCMC.packagenames="ergm",
                         # Number of proposals within each time step.
                         MCMC.burnin=100,

                         # The reason MCMC.interval=MCMC.burnin is
                         # that both represent the number of MH
                         # proposals between approximately independent
                         # draws.
                         CMLE.control.form=control.ergm(init=init.form, MCMC.prop.weights=MCMC.prop.weights.form, MCMC.prop.args=MCMC.prop.args.form, MCMC.init.maxedges=MCMC.init.maxedges, MCMC.packagenames=MCMC.packagenames, MCMC.interval=MCMC.burnin),
                         CMLE.control.diss=control.ergm(init=init.diss, MCMC.prop.weights=MCMC.prop.weights.diss, MCMC.prop.args=MCMC.prop.args.diss, MCMC.init.maxedges=MCMC.init.maxedges, MCMC.packagenames=MCMC.packagenames, MCMC.interval=MCMC.burnin),

                         EGMoME.main.method=c("Robbins-Monro"),
                         
                         SAN.maxit=10,
                         SAN.control=control.san(coef=init.form,
                           SAN.prop.weights=MCMC.prop.weights.form,
                           SAN.prop.args=MCMC.prop.args.form,
                           SAN.init.maxedges=MCMC.init.maxedges,
                           
                           SAN.burnin=MCMC.burnin,
                           SAN.packagenames=MCMC.packagenames,
                           
                           parallel=parallel,
                           parallel.type=parallel.type,
                           parallel.version.check=parallel.version.check),

                         RM.burnin=100,

                         # Plot the progress of the optimization.
                         RM.plot.progress=FALSE,
                         
                         # Initial gain --- if the process initially goes
                         # crazy beyond recovery, lower this.
                         RM.init.gain=0.1,                         
                         
                         RM.runlength=25, # Number of jumps per .C call.

                         # Interval --- number of steps between
                         # successive jumps --- is computed
                         # adaptively.
                         RM.target.ac=0.5, # Target serial autocorrelation.
                         RM.init.interval=500, # Strting interval.
                         RM.min.interval=20, # The lowest it can go.

                        
                         RM.phase1.tries=20, # Number of iterations of trying to find a reasonable configuration. FIXME: nothing happens if it's exceeded.
                         RM.phase1.jitter=0.1, # Initial jitter sd of each parameter..
                         RM.phase1.min.nonextreme=0.5, # Fraction of realizations of a statistic that are not at an extrme before it's considered "unstuck".
                         RM.phase1.max.p=0.01, # P-value that a gradient estimate must obtain before it's accepted (since sign is what's important).

                         RM.phase2sub=40, # Number of gain levels to go through.
                         RM.phase2regain=100, # Number of times gain a subphase can be repeated if the optimization is "going somewhere".
                         RM.stepdown.subphases=10, # Number of subphases to use to see whether the optimization is going somewhere.
                         RM.stepdown.p=0.5, # If the combined p-value for the trend in the parameters is less than this, repeat the subphase.
                         RM.gain.decay=0.9, # Gain decay factor.
                         RM.keep.oh=0.5, # Fraction of optimization history that is used for gradient and covariance calculation.
                         RM.jitter.mul=0.2, # The jitter standard deviation of each parameter is this times its standard deviation sans jitter.
                         RM.phase2.refine=TRUE, # Whether to use linear interpolation to refine the estimate after every run.

                         

                         RM.refine=c("linear","mean","none"), # Method, if any, used to refine the point estimate: linear interpolation, average, and none for the last value.
                         
                         RM.se=FALSE, # Whether to run Phase 3 to compute the standard errors.
                         RM.phase3n=1000, # This times the interval is the number of steps to estimate the standard errors.

                         seed=NULL,
                         parallel=0,
                         parallel.type=NULL,
                         parallel.version.check=TRUE){
  
  match.arg.pars=c("EGMoME.main.method","RM.refine")
  
  control<-list()
  formal.args<-formals(sys.function())
  for(arg in names(formal.args))
    control[[arg]]<-get(arg)

  for(arg in match.arg.pars)
    control[[arg]]<-match.arg(control[[arg]][1],eval(formal.args[[arg]]))
  
  control
}
