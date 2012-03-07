#  File ergm/R/control.stergm.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
###########################################################################
# The <control.stergm> function allows the ergm fitting process to be tuned
# by returning a list of several control parameters
############################################################################
control.stergm<-function(init.form=NULL,
                         init.diss=NULL,
                         init.method=NULL,

                         MCMC.prop.weights.form="default",MCMC.prop.args.form=NULL,
                         MCMC.prop.weights.diss="default",MCMC.prop.args.diss=NULL,
                         MCMC.init.maxedges=20000,
                         MCMC.init.maxchanges=20000,
                         MCMC.packagenames="ergm",
                         # Number of proposals within each time step.
                         MCMC.burnin=1000,

                         # The reason MCMC.interval=MCMC.burnin is
                         # that both represent the number of MH
                         # proposals between approximately independent
                         # draws.
                         CMLE.control=NULL,
                         CMLE.control.form=control.ergm(init=init.form, MCMC.prop.weights=MCMC.prop.weights.form, MCMC.prop.args=MCMC.prop.args.form, MCMC.init.maxedges=MCMC.init.maxedges, MCMC.packagenames=MCMC.packagenames, MCMC.interval=MCMC.burnin, parallel=parallel, parallel.type=parallel.type, parallel.version.check=parallel.version.check),
                         CMLE.control.diss=control.ergm(init=init.diss, MCMC.prop.weights=MCMC.prop.weights.diss, MCMC.prop.args=MCMC.prop.args.diss, MCMC.init.maxedges=MCMC.init.maxedges, MCMC.packagenames=MCMC.packagenames, MCMC.interval=MCMC.burnin, parallel=parallel, parallel.type=parallel.type, parallel.version.check=parallel.version.check),

                         EGMME.main.method=c("Stochastic-Approximation"),
                         
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

                         SA.burnin=1000,

                         # Plot the progress of the optimization.
                         SA.plot.progress=FALSE,
                         SA.max.plot.points=400,
                         
                         # Initial gain --- if the process initially goes
                         # crazy beyond recovery, lower this.
                         SA.init.gain=0.1,
                         SA.gain.decay=0.8, # Gain decay factor.
                         SA.gain.max.dist.boost=8, # Gain is boosted by the square root of the average mahalanobis distance between observed and simulated, up to this much.
                         
                         SA.runlength=25, # Number of jumps per .C call.

                         # Interval --- number of steps between
                         # successive jumps --- is computed
                         # adaptively.
                         SA.interval.mul=2, # Set the mean duration of extant ties this to be the interval.
                         SA.init.interval=500, # Starting interval.
                         SA.min.interval=20, # The lowest it can go.

                        
                         SA.phase1.tries=20, # Number of iterations of trying to find a reasonable configuration. FIXME: nothing happens if it's exceeded.
                         SA.phase1.jitter=0.1, # Initial jitter sd of each parameter..
                         SA.phase1.max.p=0.001, # P-value that a gradient estimate must obtain before it's accepted (since sign is what's important).
                         SA.phase1.backoff.rat=1.05, # If a run produces this relative increase in the objective function, it will be backed off.                         
                         SA.phase2.levels=10, # Number of gain levels to go through.
                         SA.phase2.repeats=400, # Maximum number of times gain a subphase can be repeated if the optimization is "going somewhere".
                         SA.stepdown.maxn=100, # Thin the draws for trend detection to get this many.
                         SA.stepdown.p=0.05, # If the combined p-value for the trend in the parameters is less than this, reset the subphase counter.
                         SA.stepdown.ct.base=5, # Baseline number of times in a row the p-value must be above SA.stepdown.p to reduce gain.
                         SA.stepdown.ct.subphase=1, # Number of times added to the baseline per subphase.
                         SA.phase2.backoff.rat=1.1, # If a run produces this relative increase in the objective function, it will be backed off.
                         SA.keep.oh=0.5, # Fraction of optimization history that is used for gradient and covariance calculation.
                         SA.phase2.jitter.mul=0.2, # The jitter standard deviation of each parameter is this times its standard deviation sans jitter.
                         SA.phase2.maxreljump=4, # Maximum jump per run, relative to the magnitude of other jumps in the history.
                         SA.phase2.refine=FALSE, # Whether to use linear interpolation to refine the estimate after every run. More trouble than it's worth.
                         

                         

                         SA.refine=c("linear","mean","none"), # Method, if any, used to refine the point estimate: linear interpolation, average, and none for the last value.
                         
                         SA.se=TRUE, # Whether to run Phase 3 to compute the standard errors.
                         SA.phase3.samplesize=1000, # This times the interval is the number of steps to estimate the standard errors.

                         seed=NULL,
                         parallel=0,
                         parallel.type=NULL,
                         parallel.version.check=TRUE){
  
  match.arg.pars=c("EGMME.main.method","SA.refine")

  if(!is.null(CMLE.control)) CMLE.control.form <- CMLE.control.diss <- CMLE.control
  
  control<-list()
  formal.args<-formals(sys.function())
  for(arg in names(formal.args))
    control[arg]<-list(get(arg))

  for(arg in match.arg.pars)
    control[arg]<-list(match.arg(control[[arg]][1],eval(formal.args[[arg]])))
  
  control
}
