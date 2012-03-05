\name{control.stergm}
\alias{control.stergm}
\title{Auxiliary for Controlling Separable Temporal ERGM Fitting
}
\description{  Auxiliary function as user interface for fine-tuning 'stergm' fitting.

}
\usage{
 control.stergm (init.form=NULL,
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

                         EGMME.main.method=c("Robbins-Monro"),
                         
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

                         RM.burnin=1000,

                         # Plot the progress of the optimization.
                         RM.plot.progress=FALSE,
                         RM.max.plot.points=400,
                         
                         # Initial gain --- if the process initially goes
                         # crazy beyond recovery, lower this.
                         RM.init.gain=0.01,                         
                         
                         RM.runlength=25, # Number of jumps per .C call.

                         # Interval --- number of steps between
                         # successive jumps --- is computed
                         # adaptively.
                         RM.target.ac=0.5, # Target serial autocorrelation.
                         RM.init.interval=500, # Starting interval.
                         RM.min.interval=20, # The lowest it can go.

                        
                         RM.phase1.tries=20, # Number of iterations of trying to find a reasonable configuration. FIXME: nothing happens if it's exceeded.
                         RM.phase1.jitter=0.1, # Initial jitter sd of each parameter..
                         RM.phase1.max.p=0.01, # P-value that a gradient estimate must obtain before it's accepted (since sign is what's important).
                         RM.phase1.backoff.rat=1.05, # If a run produces this relative increase in the objective function, it will be backed off.
                         
                         RM.phase2.levels=20, # Number of gain levels to go through.
                         RM.phase2.repeats=400, # Maximum number of times gain a subphase can be repeated if the optimization is "going somewhere".
                         RM.stepdown.runs=40, # Number of subphases to use to see whether the optimization is going somewhere.
                         RM.stepdown.thin=10, # Thin the draws for trend detection.
                         RM.stepdown.p=0.01, # If the combined p-value for the trend in the parameters is less than this, reset the subphase counter.
                         RM.stepdown.ct=5, # Number of times in a row the p-value must be above RM.stepdown.p to reduce gain.
                         RM.gain.decay=0.7, # Gain decay factor.
                         RM.keep.oh=0.75, # Fraction of optimization history that is used for gradient and covariance calculation.
                         RM.phase2.jitter.mul=0.2, # The jitter standard deviation of each parameter is this times its standard deviation sans jitter.
                         RM.phase2.maxreljump=16, # Maximum jump per run, relative to the magnitude of other jumps in the history.
                         RM.phase2.refine=TRUE, # Whether to use linear interpolation to refine the estimate after every run.
                         

                         

                         RM.refine=c("linear","mean","none"), # Method, if any, used to refine the point estimate: linear interpolation, average, and none for the last value.
                         
                         RM.se=TRUE, # Whether to run Phase 3 to compute the standard errors.
                         RM.phase3.samplesize=1000, # This times the interval is the number of steps to estimate the standard errors.

                         seed=NULL,
                         parallel=0,
                         parallel.type=NULL,
                         parallel.version.check=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{init.form, init.diss}{numeric or \code{NA} vector equal in length to the number of
    parameters in the formation/dissolution model or \code{NULL} (the default); the initial
    values for the estimation and coefficient offset terms. If \code{NULL}
    is passed, all of the initial values are computed using the method specified by
    \code{\link[=control.ergm]{control$init.method}}.
    If a numeric vector is given, the elements of the vector are interpreted as follows:
    \itemize{
      \item{Elements corresponding to
        terms enclosed in \code{offset()} are used as the fixed offset
        coefficients.
      }
      \item{Elements that do not correspond to offset terms and are not
        \code{NA} are used as starting values in the estimation. 
      }
      \item{Initial values for the elements that are \code{NA} are fit
        using the method specified by
        \code{\link[=control.ergm]{control$init.method}}.
      }
    }
    Passing \code{control.ergm(init=coef(prev.fit))} can be used to ``resume''
    an uncoverged \code{\link{ergm}} run, but see
    \code{\link{enformulate.curved}}.
  }
  \item{init.method}{Estimation method
    used to acquire initial values for estimation. Unused at this time.
    }
  \item{MCMC.prop.weights.form, MCMC.prop.weights.diss}{Specifies the method to allocate probabilities of
    being proposed to dyads in the formation/dissolution phase. Defaults to \code{"default"}, which picks a
    reasonable default for the specified constraint. Possible values include
    \code{"TNT"}, \code{"random"}, though not
    all values may be used
    with all possible constraints.
  }
  \item{MCMC.prop.args.form, MCMC.prop.args.diss}{
    An alternative, direct way of specifying additional arguments to
    the proposal in the formation/dissolution phase.
  }
  \item{MCMC.init.maxedges}{Maximum number 
    of edges for which to allocate space.
  }
  \item{MCMC.init.maxchanges}{Maximum number of changes in dynamic network
    simulation for which to allocate space.
  }
  \item{MCMC.packagenames}{Names of packages in which change statistic
    functions are found.
  }
  \item{MCMC.burnin}{Number of Metropolis-Hastings steps per phase
                        (formation and dissolution) per time step.
  }
  \item{CMLE.control}{A convenience argument for specifying both
    \code{CMLE.control.form} and \code{CMLE.control.diss} at once.  See
    \code{\link{control.ergm}}.
  }
  \item{CMLE.control.form, CMLE.control.diss}{Control parameters used to
                         fit the CMLE for the formation/dissolution ERGM.  See
    \code{\link{control.ergm}}.
  }
  \item{EGMME.main.method}{Estimation method
    used to find the Equilibrium Generalized Method of Moments
    estimator.  Currently only "Robbins-Monro" is implemented.
  }
  \item{SAN.maxit}{ When \code{target.stats} argument is passed to
  \code{\link{ergm}}, the maximum number of attempts to use
  \code{\link{san}} to obtain a network with statistics close to those
  specified.
  }
  \item{SAN.control}{SAN control parameters.  See
    \code{\link{control.san}}
  }
  \item{RM.burnin}{Number of time steps to advance the starting network
    before beginning the optimization.
  }
  \item{RM.plot.progress}{Logical: If TRUE, plot the trajectories of the
                         parameters and target statistics as the
                         optimization progresses. Do NOT use with non-interactive plotting devices
    like \code{\link{pdf}}.
  }

  \item{RM.max.plot.points}{
    If \code{RM.plot.progress==TRUE}, the maximum number of time points
    to be plotted. Defaults to 400. If more iterations elapse, they will
    be thinned to at most 400 before plotting.
    }
  \item{RM.init.gain}{Initial gain, the multiplier for the parameter update size.  If the process initially goes
    crazy beyond recovery, lower this value.
  }
  \item{RM.runlength}{Number of parameter trials and updates per C run.
  }
  \item{RM.target.ac}{Target serial autocorrelation. The number of time
    steps between updates of the parameters will be increased and
    decreased to target this value.
  }
  \item{RM.init.interval}{Initial number of time steps between updates
    of the parameters.
  }
  \item{RM.min.interval}{The lowest the number of time steps between updates
    of the parameters can go.
  }
  
  \item{RM.phase1.tries}{Number of runs trying to find a 
    reasonable parameter and network configuration.
  }
  \item{RM.phase1.jitter}{Initial jitter standard deviation
    of each parameter.
  }
  \item{RM.phase1.max.p}{P-value that a gradient estimate must 
    obtain before it is accepted (since sign is what is important).
  }

  \item{RM.phase1.backoff.rat}{If the run produces this 
    relative increase in the approximate objective function, it will be backed off.
  }
  \item{RM.phase2.levels}{Number of gain levels (subphases) to go through.
  }
  \item{RM.phase2.repeats, RM.stepdown.runs, RM.stepdown.thin,
    RM.stepdown.p, RM.stepdown.ct}{A gain level may be repeated multiple times (up to
    \code{RM.phase2.repeats}) if the optimizer detects that the
    parameters are "going somewhere": moving in a clear direction, so
    slowing them down at that point is counterproductive. To detect this
    trend, it looks at parameter values up to \code{RM.phase2.stepdown.runs} back, thinning them
    by \code{RM.stepdown.thin}, and fitting a GLS model for a linear trend, with AR(2)
    autocorrelation for each variable, combining the p-values. If this
    combined p-value for the trend in at least one parameter is greater
    than \code{RM.stepdown.p} \code{RM.stepdown.ct} runs in a row, the
    gain level (subphase) is allowed to end. Otherwise, the process
    continues at the same gain level.
  }
  \item{RM.gain.decay}{Gain decay factor.
  }
  \item{RM.keep.oh}{Fraction of optimization history
    used for gradient and covariance estimation.
  }
  \item{RM.phase2.jitter.mul}{Jitter standard deviation of each 
    parameter is this value times its standard deviation without jitter.
  }
  \item{RM.phase2.refine}{Logical: If TRUE, use linear 
    interpolation to refine the estimate after every run.
  }

  \item{RM.phase2.maxreljump}{To keep the optimization from "running
  away" due to, say, a poor gradient estimate building on itself, if a
  magnitude of change (Mahalanobis distance) in parameters over the course of a run divided by
  average magnitude of change for recent runs exceeds this, the change
  is truncated to this amount times the average for recent runs.}
  

  \item{RM.refine}{Method, if any, used to 
    refine the point estimate at the end: "linear" for linear interpolation, 
    "mean" for average, and "none" to use the last value.
  }
  \item{RM.se}{Logical:  If TRUE (the default), get an MCMC sample of statistics at
    the final estimate and compute the 
    covariance matrix (and hence standard errors) of the
    parameters. This sample is stored and can also be used by
    \code{\link{mcmc.diagnostics.stergm}} to assess convergence.
  }
  \item{RM.phase3.samplesize}{This multiplied by the interval is the 
    number of steps to estimate the standard errors.
  }
  \item{seed}{Seed value (integer) for the random number generator.  
    See \code{\link[base]{set.seed}}
  }
  \item{parallel}{Number of threads in which to run the
    sampling. Defaults to 0 (no parallelism). See the entry on
    \link[=ergm-parallel]{parallel processing} for
    details and troubleshooting.
  }
  \item{parallel.type}{API to use for parallel
    processing. Supported values are \code{"MPI"} and
    \code{"SOCK"}. Defaults to using the \code{snow} package
    default.
  }
  \item{parallel.version.check}{Logical:  If TRUE, check that the version of
    \code{\link[=ergm-package]{ergm}} running on the slave nodes is the
    same as that running on the master node.
  }
}
\details{
      This function is only used within a call to the \code{\link{stergm}} function.
    See the \code{usage} section in \code{\link{stergm}} for details.
}
\value{
  A list with arguments as components.
}
\references{
    \itemize{
    \item
    Boer, P., Huisman, M., Snijders, T.A.B., and Zeggelink, E.P.H. (2003),
    StOCNET User\'s Manual. Version 1.4.
    
    \item Firth (1993),
    Bias Reduction in Maximum Likelihood Estimates.
    Biometrika,
    80: 27-38.

    \item Hunter, D. R. and M. S. Handcock (2006), Inference in curved
    exponential family models for networks. Journal of Computational
    and Graphical Statistics, 15: 565-583.
    
    \item Hummel, R. M., Hunter, D. R., and Handcock, M. S. (2010),
    A Steplength Algorithm for Fitting ERGMs, Penn State Department
    of Statistics Technical Report.
  }
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{stergm}}. The \code{\link{control.simulate.stergm}} 
function performs a 
similar function for
\code{\link{simulate.stergm}}.  
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
