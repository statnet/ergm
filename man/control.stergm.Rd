%  File ergm/man/control.stergm.Rd
%  Part of the statnet package, http://statnet.org
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
%  Copyright 2012 the statnet development team
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{control.stergm}
\alias{control.stergm}
\title{Auxiliary for Controlling Separable Temporal ERGM Fitting
}
\description{  Auxiliary function as user interface for fine-tuning 'stergm' fitting.

}
\usage{
          control.stergm(init.form=NULL,
                         init.diss=NULL,
                         init.method=NULL,
                         MCMC.prop.weights.form="default",MCMC.prop.args.form=NULL,
                         MCMC.prop.weights.diss="default",MCMC.prop.args.diss=NULL,
                         MCMC.init.maxedges=20000,
                         MCMC.init.maxchanges=20000,
                         MCMC.packagenames="ergm",
                         MCMC.burnin=1000,
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
                         SA.plot.progress=FALSE,
                         SA.max.plot.points=400,
                         SA.init.gain=0.1,
                         SA.gain.decay=0.8,
                         SA.gain.max.dist.boost=8,
                         SA.runlength=25,
                         SA.interval.mul=2,
                         SA.init.interval=500, 
                         SA.min.interval=20, 
                         SA.phase1.tries=20, 
                         SA.phase1.jitter=0.1, 
                         SA.phase1.max.p=0.001, 
                         SA.phase1.backoff.rat=1.05, 
                         SA.phase2.levels=10, 
                         SA.phase2.repeats=400,
                         SA.stepdown.maxn=100,
                         SA.stepdown.p=0.05,
                         SA.stepdown.ct.base=5,
                         SA.stepdown.ct.subphase=1,
                         SA.phase2.backoff.rat=1.1, 
                         SA.keep.oh=0.5, 
                         SA.phase2.jitter.mul=0.2, 
                         SA.phase2.maxreljump=4, 
                         SA.phase2.refine=FALSE, 
                         SA.refine=c("linear","mean","none"), 
                         SA.se=TRUE, 
                         SA.phase3.samplesize=1000, 
                         seed=NULL,
                         parallel=0,
                         parallel.type=NULL,
                         parallel.version.check=TRUE)
}
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
  \item{SA.burnin}{Number of time steps to advance the starting network
    before beginning the optimization.
  }
  \item{SA.plot.progress}{Logical: If TRUE, plot the trajectories of the
                         parameters and target statistics as the
                         optimization progresses. Do NOT use with non-interactive plotting devices
    like \code{\link{pdf}}.
  }

  \item{SA.max.plot.points}{
    If \code{SA.plot.progress==TRUE}, the maximum number of time points
    to be plotted. Defaults to 400. If more iterations elapse, they will
    be thinned to at most 400 before plotting.
    }
  \item{SA.init.gain}{Initial gain, the multiplier for the parameter update size.  If the process initially goes
    crazy beyond recovery, lower this value.
  }
  \item{SA.gain.decay}{Gain decay factor.
  }
  \item{SA.gain.max.dist.boost}{
    Gain is boosted or shrunk by the square root of the average mahalanobis
    distance between observed and simulated, up to this much. This
    speeds up convergence once the process finds the right values.
    }
  \item{SA.runlength}{Number of parameter trials and updates per C run.
  }
  \item{SA.interval.mul}{The number of time steps between updates of the
    parameters is set to be this times the mean duration of extant ties.}
  \item{SA.init.interval}{Initial number of time steps between updates
    of the parameters.
  }
  \item{SA.min.interval}{The lowest the number of time steps between updates
    of the parameters can go.
  }
  
  \item{SA.phase1.tries}{Number of runs trying to find a 
    reasonable parameter and network configuration.
  }
  \item{SA.phase1.jitter}{Initial jitter standard deviation
    of each parameter.
  }
  \item{SA.phase1.max.p}{P-value that a gradient estimate must 
    obtain before it is accepted (since sign is what is important).
  }

  \item{SA.phase1.backoff.rat, SA.phase2.backoff.rat}{If the run produces this 
    relative increase in the approximate objective function, it will be backed off.
  }
  \item{SA.phase2.levels}{Number of gain levels (subphases) to go through.
  }
  \item{SA.phase2.repeats, SA.stepdown.maxn,
    SA.stepdown.p, SA.stepdown.ct.base, SA.stepdown.ct.subphase}{A gain level may be repeated multiple times (up to
    \code{SA.phase2.repeats}) if the optimizer detects that the
    objective function is improving, so slowing down the parameters at that point is counterproductive. To detect this
    it looks at the the window controlled by \code{SA.keep.oh}, thinning
    objective function values
    to get \code{SA.stepdown.maxn}, and fitting a GLS model for a linear trend, with AR(2)
    autocorrelation. If there is a negative trend is greater than
    \code{SA.stepdown.p} \code{SA.stepdown.ct.base +
    SA.stepdown.ct.subphase * subphase} runs in a row, the
    gain level (subphase) is allowed to end. Otherwise, the process
    continues at the same gain level.
  }
  \item{SA.keep.oh}{Fraction of optimization history
    used for gradient and covariance estimation.
  }
  \item{SA.phase2.jitter.mul}{Jitter standard deviation of each 
    parameter is this value times its standard deviation without jitter.
  }
  \item{SA.phase2.refine}{Logical: If TRUE, use linear 
    interpolation to refine the estimate after every run.
  }

  \item{SA.phase2.maxreljump}{To keep the optimization from "running
  away" due to, say, a poor gradient estimate building on itself, if a
  magnitude of change (Mahalanobis distance) in parameters over the course of a run divided by
  average magnitude of change for recent runs exceeds this, the change
  is truncated to this amount times the average for recent runs.}
  

  \item{SA.refine}{Method, if any, used to 
    refine the point estimate at the end: "linear" for linear interpolation, 
    "mean" for average, and "none" to use the last value.
  }
  \item{SA.se}{Logical:  If TRUE (the default), get an MCMC sample of statistics at
    the final estimate and compute the 
    covariance matrix (and hence standard errors) of the
    parameters. This sample is stored and can also be used by
    \code{\link{mcmc.diagnostics.stergm}} to assess convergence.
  }
  \item{SA.phase3.samplesize}{This multiplied by the interval is the 
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


\seealso{
  \code{\link{stergm}}. The \code{\link{control.simulate.stergm}} 
function performs a 
similar function for
\code{\link{simulate.stergm}}.  
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
