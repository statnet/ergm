\name{control.stergm}
\alias{control.stergm}
\title{Auxiliary for Controlling Separable Temporal ERGM Fitting
}
\description{  Auxiliary function as user interface for fine-tuning 'stergm' fitting.

}
\usage{
control.stergm<-function(init.form=NULL,
                         init.diss=NULL,
                         init.method=NULL,

                         MCMC.prop.weights.form="default",
                         MCMC.prop.args.form=NULL,
                         MCMC.prop.weights.diss="default",
                         MCMC.prop.args.diss=NULL,
                         MCMC.init.maxedges=20000,
                         MCMC.init.maxchanges=20000,
                         MCMC.packagenames="ergm",
                         MCMC.burnin=1000,

                         CMLE.control.form=control.ergm(init=init.form,
                                             MCMC.prop.weights=MCMC.prop.weights.form, 
                                             MCMC.prop.args=MCMC.prop.args.form, 
                                             MCMC.init.maxedges=MCMC.init.maxedges, 
                                             MCMC.packagenames=MCMC.packagenames, 
                                             MCMC.interval=MCMC.burnin),
                         CMLE.control.diss=control.ergm(init=init.diss, 
                                             MCMC.prop.weights=MCMC.prop.weights.diss, 
                                             MCMC.prop.args=MCMC.prop.args.diss, 
                                             MCMC.init.maxedges=MCMC.init.maxedges, 
                                             MCMC.packagenames=MCMC.packagenames, 
                                             MCMC.interval=MCMC.burnin),

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

                         RM.burnin=1000,
                         RM.plot.progress=FALSE,
                         RM.init.gain=0.1,                         
                         RM.runlength=25, 
                         RM.target.ac=0.5, 
                         RM.init.interval=500,
                         RM.min.interval=20,
                         RM.phase1.tries=20,
                         RM.phase1.jitter=0.1,
                         RM.phase1.min.nonextreme=0.5,
                         RM.phase1.max.p=0.01,
                         RM.phase2sub=20,
                         RM.phase2regain=100,
                         RM.stepdown.subphases=10,
                         RM.stepdown.p=0.5,
                         RM.gain.decay=0.7,
                         RM.keep.oh=0.5,
                         RM.jitter.mul=0.2,
                         RM.phase2.refine=TRUE,
                         RM.phase2.refine.maxjump=3,
                         RM.refine=c("linear","mean","none"),
                         RM.se=FALSE,
                         RM.phase3n=1000,

                         seed=NULL,
                         parallel=0,
                         parallel.type=NULL,
                         parallel.version.check=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  
  \item{init.form}{numeric or \code{NA} vector equal in length to the number of
    parameters in the formation model or \code{NULL} (the default); the initial
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
  \item{init.diss}{As \code{init.diss} for the dissolution model.
  }
  \item{init.method}{Unused for now.  The default method for finding the
    starting coefficient values (if \code{init} is not specified) is 
    maximum pseudo-likelihood estimation (MPLE).
  }
  \item{MCMC.prop.weights.form}{Specifies the method to allocate probabilities of
    being proposed to dyads in the formation phase. Defaults to \code{"default"}, which picks a
    reasonable default for the specified constraint. Possible values are
    \code{"TNT"}, \code{"random"}, and \code{"nonobserved"}, though not
    all values may be used
    with all possible constraints (in the \code{\link{ergm}} function).
  }
  \item{MCMC.prop.args.form}{
    An alternative, direct way of specifying additional arguments to
    the proposal in the formation phase.
  }
  \item{MCMC.prop.weights.diss}{As \code{MCMC.prop.weights.form}, 
    for the dissolution phase.
  }
  \item{MCMC.prop.args.diss}{As \code{MCMC.prop.args.form}, 
    for the dissolution phase.
  }
  \item{MCMC.maxedges}{Maximum number 
    of edges for which to allocate space.
  }
  \item{MCMC.maxchanges}{Maximum number of changes in dynamic network
    simulation for which to allocate space.
  }
  \item{MCMC.packagenames="ergm",}{ Not yet documented. 
  }
  \item{MCMC.burnin}{Number of proposals within each time step.
    This number functions as both \code{MCMC.burnin} and
    \code{MCMC.interval}.  They are the same because each
    represents the number of MH proposals between approximately
    independent draws.
  }
  \item{CMLE.control.form}{Control parameters used to fit 
    cross-sectional ERGMs using the formation model.  See
    \code{\link{control.ergm}}.
  }
  \item{CMLE.control.diss}{Control parameters used to fit 
    cross-sectional ERGMs using the dissolution model.  See
    \code{\link{control.ergm}}.
  }
  \item{EGMoME.main.method}{Estimation method
    used to find the Equilibrium Generalized Method of Moments
    estimator.  Currently only "Robbins-Monro" is implemented.
  }
  \item{SAN.maxit}{ Not yet documented. 
  }
  \item{SAN.control}{SAN control parameters.  See
    \code{\link{control.san}}
  }
  \item{RM.burnin}{Not yet documented. 
  }
  \item{RM.plot.progress}{Logical: If TRUE, plot the progress of
    the optimization.
  }
  \item{RM.init.gain}{Initial gain.  If the process initially goes
    crazy beyond recovery, lower this value.
  }
  \item{RM.runlength}{Number of jumps per .C call.
  }
  \item{RM.target.ac}{Target serial autocorrelation.
  }
  \item{RM.init.interval}{Strating interval.  The interval,
    or number of steps between successive jumps,
    is computed adaptively.
  }
  \item{RM.min.interval}{The lowest the interval can go.
  }
  \item{RM.phase1.tries}{Number of iterations of trying to find a 
    reasonable configuration. (Currently not used;
    nothing happens if it's exceeded.)
  }
  \item{RM.phase1.jitter}{Initial jitter standard deviation
    of each parameter.
  }
  \item{RM.phase1.min.nonextreme}{Fraction of realizations of a 
    statistic that are not at an extrme before it's considered "unstuck".
  }
  \item{RM.phase1.max.p}{P-value that a gradient estimate must 
    obtain before it is accepted (since sign is what is important).
  }
  \item{RM.phase2sub}{Number of gain levels to go through.
  }
  \item{RM.phase2regain}{Number of times gain a subphase can be 
    repeated if the optimization is "going somewhere".
  }
  \item{RM.stepdown.subphases}{Number of subphases to use to see 
    whether the optimization is going somewhere.
  }
  \item{RM.stepdown.p}{If the combined p-value for the trend in 
    the parameters is less than this value, repeat the subphase.
  }
  \item{RM.gain.decay}{Gain decay factor.
  }
  \item{RM.keep.oh}{Fraction of optimization history
    used for gradient and covariance calculation.
  }
  \item{RM.jitter.mul}{Jitter standard deviation of each 
    parameter is this value times its standard deviation without jitter.
  }
  \item{RM.phase2.refine}{Logical: If TRUE, use linear 
    interpolation to refine the estimate after every run.
  }
  \item{RM.phase2.refine.maxjump}{Maximum linear interpolator jump 
    magnitude (in standard deviations of the parameter from the 
    current mean).
  }
  \item{RM.refine=c("linear","mean","none"), }{Method, if any, used to 
    refine the point estimate: "linear" for linear interpolation, 
    "mean" for average, and "none" for the last value.
  }
  \item{RM.se}{Logical:  If TRUE, run Phase 3 to compute the 
    standard errors.
  }
  \item{RM.phase3n}{This multiplied by the interval is the 
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
  
  

% Below is some old documentation that may be helpful to augment what is
% currently given:
%
%  \item{style}{character; The style of method of moments estimation
%    to use. The default is a form of stochastic approximation
%    (\code{"Robbins-Monro"}), but it should only be used if it is known
%    a priori that the derivative of each element of the equilibrium
%    expected values of statistics of interest with respect to the
%    corresponding formation phase parameter is positive. The other
%    option, (\code{"SPSA"}) is less precise but does not make this
%    assumption. \code{"SPSA2"} is multithreaded, but has more overhead,
%    and may give wrong results, since \R's routines are not threadsafe.
%}
%  \item{RM.phase1n_base}{this helps define the 'phase1n' param, which in turn
%                     multiplies 'RM.interval' to control the number of
%                     phase1 iterations; this is the base portion of 'phase1n',
%                     which is added to 3*(the number of formation coefficients)
%}
%  \item{RM.phase2n_base}{phase2 is a 3-deep nested for-loop and 'RM.phase2sub' limits
%}
%\item{RM.phase2sub}{
%  phase2 is a 3-deep nested for-loop and 'RM.phase2sub' limits
%}
%\item{RM.init_gain}{
%  this is only used to adjust 'aDdiaginv'in phase1,
%                      in particular:
%                             aDdiaginv = gain/sqrt(aDdiaginv)
%                      default=0.5
%}
%
%  \item{RM.interval}{like the SPSA.interval, this seems a little more like
%                      a sample size, than an interval, it helps control the
%                      number of MCMCsteps used in phase1 and phase2; in
%                      phase2, this limits the innermost loop counter; default=100
%}
%  \item{RM.burnin}{the number of MCMC steps to disregard for the burn-in
%                      period; default=1000
%}
%  \item{SPSA.a,SPSA.alpha,SPSA.A,SPSA.c,SPSA.gamma}{ Controls the rate
%                      of decay of SPSA gain as
%                           SPSA.c/(i+1)^(SPSA.gamma)
%                      where i is indexed from 0 to SPSA.iterations
%}
%  \item{SPSA.iterations}{the number of iterations to use in the SPSA optimization;
%                      default=1000
%}
%  \item{SPSA.interval}{this is eventually received as 'S' and looks like a
%                      a sample size, rather than an interval, since 'S' controls
%                      the number of MCMC steps that contribute to the stats vector
%}
%  \item{SPSA.burnin}{the number of MCMC steps to disregard for the burnin
%                      period; default=1000
%}
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
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
