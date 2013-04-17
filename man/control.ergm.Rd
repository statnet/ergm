%  File man/control.ergm.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2013 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{control.ergm}
\alias{control.ergm}
\title{ Auxiliary for Controlling ERGM Fitting }
\description{
  Auxiliary function as user interface for fine-tuning 'ergm' fitting.
}
\usage{
control.ergm(drop=TRUE,

             init=NULL,
             init.method=NULL,
                       
             main.method=c("MCMLE","Robbins-Monro",
               "Stochastic-Approximation","Stepping"),
             force.main=FALSE,
             main.hessian=TRUE,

             MPLE.max.dyad.types=1e+6, 
             MPLE.samplesize=50000,                       
             MPLE.type=c("glm", "penalized"),
                      
             MCMC.prop.weights="default", 
             MCMC.prop.args=list(),
             MCMC.burnin=10000,
             MCMC.interval=100,
             MCMC.samplesize=10000,
             MCMC.return.stats=TRUE,
             MCMC.burnin.retries=0,
             MCMC.burnin.check.last=1/2,
             MCMC.burnin.check.alpha=0.01,
             MCMC.runtime.traceplot=FALSE,
             MCMC.init.maxedges=20000,
             MCMC.max.maxedges=Inf,
             MCMC.addto.se=TRUE,
             MCMC.compress=FALSE,
             MCMC.packagenames=c(),

             SAN.maxit=10,
             SAN.control=control.san(
                         coef=init,
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
             obs.MCMC.samplesize=MCMC.samplesize,
             obs.MCMC.interval=MCMC.interval,
             obs.MCMC.burnin=MCMC.burnin,
             MCMLE.check.degeneracy=FALSE,
             MCMLE.MCMC.precision=0.05,
             MCMLE.metric=c("lognormal", "Median.Likelihood",
               "EF.Likelihood", "naive"),
             MCMLE.method=c("BFGS","Nelder-Mead"),
             MCMLE.trustregion=20,
             MCMLE.dampening=FALSE,
             MCMLE.dampening.min.ess=20,
             MCMLE.dampening.level=0.1,
             MCMLE.steplength=0.5,
             MCMLE.adaptive.trustregion=3,
             MCMLE.adaptive.epsilon=0.01,
             MCMLE.sequential=TRUE,
             MCMLE.density.guard.min=10000,
             MCMLE.density.guard=exp(3),


             SA.phase1_n=NULL, 
             SA.initial_gain=NULL, 
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

             loglik.control=control.logLik.ergm(),

             seed=NULL,
             parallel=0,
             parallel.type=NULL,
             parallel.version.check=TRUE,
             
             \dots)
}
\details{
    This function is only used within a call to the \code{\link{ergm}} function.
    See the \code{usage} section in \code{\link{ergm}} for details.
}
\arguments{
  \item{drop}{Logical:  If TRUE, terms whose observed statistic values
    are at the extremes of their possible ranges are dropped from the
    fit and their corresponding parameter estimates are set to
    plus or minus infinity, as appropriate.  This is done because 
    maximum likelihood estimates cannot exist when the vector of
    observed statistic lies on the boundary of the convex hull of
    possible statistic values.
  }
  \item{init}{numeric or \code{NA} vector equal in length to the number of
    parameters in the model or \code{NULL} (the default); the initial
    values for the estimation and coefficient offset terms. If \code{NULL}
    is passed, all of the initial values are computed using the method specified by
    \code{\link[=control.ergm]{control$init.method}}.
    If a numeric vector is given, the elements of the vector are interpreted as follows:
    \itemize{
      \item{Elements corresponding to
        terms enclosed in \code{offset()} are used as the fixed offset
        coefficients. Note that offset coefficients alone can be more conveniently specified
	using  \code{\link{ergm}} argument \code{offset.coef}. If both
	\code{offset.coef} and \code{init} arguments are given,
	values in \code{offset.coef} will take precedence.
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
  \item{init.method}{A chatacter vector or \code{NULL}. The default method for 
    finding the starting coefficient values, if \code{init} is not specified, is
    maximum pseudo-likelihood estimation (MPLE). Another valid value is
    \code{"zeros"} for a vector of \code{0} of appropriate length.
  }
  \item{main.method}{One of "MCMLE","Robbins-Monro",
    "Stochastic-Approximation", or "Stepping".  Chooses
    the estimation method used to find the MLE.
    \code{MCMLE} attempts to maximize an approximation to the
    log-likelihood function.  \code{Robbins-Monro} and
    \code{Stochastic-Approximation} are both 
    stochastic approximation algorithms that try to 
    solve the method of moments equation that yields the
    MLE in the case of an exponential family model.
    Another alternative is a partial stepping algorithm
    (\code{Stepping}) as in Hummel et al. (2012).
    The direct use of the likelihood function has many theoretical
    advantages over stochastic approximation, but the choice will
    depend on the model and data being fit. See Handcock (2000) and 
    Hunter and Handcock (2006) for details.
  }
  \item{force.main}{Logical:  If TRUE, then force MCMC-based estimation
    method, even if the exact MLE can be computed via maximum pseudolikelihood
    estimation.
  }
  \item{main.hessian}{Logical:  If TRUE, then an approximate Hessian
    matrix is used in the MCMC-based estimation method.
  }
  \item{MPLE.max.dyad.types}{Maximum number of unique values of 
    change statistic vectors, which are the predictors in a logistic
    regression used to calculate the MPLE.  This calculation uses 
    a compression algorithm that allocates space based on 
    \code{MPLE.max.dyad.types}.
  }
  \item{MPLE.samplesize}{Not currently documented; used in 
    conditional-on-degree version of MPLE.
  }
  \item{MPLE.type}{One of "glm" or "penalized".  Chooses method of
    calculating MPLE.  "glm" is the usual
    formal logistic regression, whereas "penalized" uses the 
    bias-reduced method
    of Firth (1993) as originally implemented by 
    Meinhard Ploner, Daniela Dunkler, Harry
    Southworth, and Georg Heinze in the "logistf" package.
  }
  \item{MCMC.prop.weights}{Specifies the proposal distribution used in the MCMC
    Metropolis-Hastings algorithm.  Possible choices are \code{"TNT"} or 
    \code{"random"}; the \code{"default"} is one of these two, depending on the
    constraints in place (as defined by the \code{constraints}
    argument of the \code{\link{ergm}} function), though not all weights
    may be used with all constraints.  
    The \code{TNT} (tie / no tie) option puts roughly equal weight on selecting a
    dyad with or without a tie as a candidate for toggling, whereas the
    \code{random} option puts equal weight on all possible dyads, though the
    interpretation of \code{random} may change according to the constraints in
    place.  When no constraints are in place, the default is TNT, which 
    appears to improve Markov chain mixing particularly for networks with a
    low edge density, as is typical of many realistic social networks.
  }
  \item{MCMC.prop.args}{An alternative, direct way of specifying additional arguments to proposal. 
  }
  \item{MCMC.burnin}{Number of proposals before any MCMC sampling
    is done. It typically is set to a fairly large number.
  }
  \item{MCMC.interval}{Number of proposals between sampled statistics.
  }
  \item{MCMC.samplesize}{Number of network statistics,
    randomly drawn from a given distribution on the set of all networks,
    returned by the Metropolis-Hastings algorithm.
  }
  \item{MCMC.return.stats}{Logical:  If TRUE, return the matrix
    of MCMC-sampled network statistics.  This matrix should have
    \code{MCMC.samplesize} rows.
    This matrix can be 
    used directly by the \code{coda} package to assess MCMC
    convergence.
  }
  \item{MCMC.burnin.retries}{Maximum number of times to rerun the burn-in phase if a failure to
    converge is detected. Defaults to 0 for no checks.
  }
  \item{MCMC.burnin.check.last}{What fraction at the end of the burn-in phase
    to use for detecting non-convergence. Defaults to one half.
  }
  \item{MCMC.burnin.check.alpha}{Maximum Bonferroni-adjusted P-value under the
    Geweke test for a statistic to be considered unconverged. Defaults
    to 0.01.
  }
  \item{MCMC.runtime.traceplot}{Logical:  If TRUE, plot traceplots of 
    the MCMC sample after every MCMC MLE iteration.
  }
  \item{MCMC.init.maxedges, MCMC.max.maxedges}{Maximum number of edges
    expected in network. Starting at \code{MCMC.init.maxedges}, it will
    be incremented by a factor of 10 if exceeded during fitting, up to
    \code{MCMC.max.maxedges}, at which point the process will stop with
    an error.
  }
  \item{MCMC.addto.se}{Not yet documented.
  }
  \item{MCMC.compress}{Logical:  If TRUE, the matrix of sample statistics
    returned is compressed to the set of unique statistics with a 
    column of frequencies post-pended.
  }
  \item{MCMC.packagenames}{Names of packages in which to look for change statistic
    functions in addition to those autodetected. This argument should not be
    needed outside of very strange setups.
  }

  \item{SAN.maxit}{When \code{target.stats} argument is passed to
  \code{\link{ergm}}, the maximum number of attempts to use
  \code{\link{san}} to obtain a network with statistics close to those
  specified.}  

  \item{SAN.control}{Control arguments to \code{\link{san}}.
    See \code{\link{control.san}} for details.
  }

  \item{MCMLE.maxit}{Maximum number of times the parameter for the MCMC
    should be updated by maximizing the MCMC likelihood. At each
    step the parameter is changed to the values that maximizes
    the MCMC likelihood based on the current sample.
  }

  \item{MCMLE.conv.min.pval}{ After every MCMC sample, a Hotelling's T^2
    test for equality of MCMC-simulated network statistics to observed
    is conducted, and if its P-value exceeds this setting, the
    estimation is considered to have converged and finishes. To turn
    this off and perform all
    \code{\link[=control.ergm]{control$MCMLE.maxit}} iterations, set
    \code{MCMLE.conv.min.pval=1}.  }
  
  \item{MCMLE.NR.maxit}{Maximum number of iterations in the
    Newton-Raphson optimization.
  }
  \item{MCMLE.NR.reltol}{Not yet documented.
  }
  \item{obs.MCMC.samplesize,obs.MCMC.burnin,obs.MCMC.interval}{
    Sample size, burnin, and interval parameters for the MCMC sampling used
    when unobserved data are present in the estimation routine.    
  }
  \item{MCMLE.check.degeneracy}{Logical:  If TRUE, employ
    a check for model degeneracy.
  }
  \item{MCMLE.MCMC.precision}{Vector of upper bounds on the precisions 
    of the standard errors induced by the MCMC algorithm.
  }
  \item{MCMLE.metric}{Method to calculate the loglikelihood approximation.
    See Hummel et al (2010)
    for an explanation of "lognormal" and "naive".
  }
  \item{MCMLE.method}{Name of the optimization method
    to use. See \code{\link[stats]{optim}} for the options. 
    The default method
    \code{"BFGS"} is a quasi-Newton method (also known as a variable
    metric algorithm). It is attributed to
    Broyden, Fletcher, Goldfarb and Shanno. This uses function values
    and gradients to build up a picture of the surface to be optimized.
  }
  \item{MCMLE.trustregion}{Maximum increase the algorithm will
    allow for the approximated likelihood at a given iteration.
    See Snijders (2002) for details.
  }
  \item{MCMLE.dampening}{(logical) Should likelihood dampening be used?}
  \item{MCMLE.dampening.min.ess}{The effective sample size below which dampening is used.}
  \item{MCMLE.dampening.level}{The proportional distance from boundary of the convex hull move.}
  \item{MCMLE.steplength}{Multiplier for step length, which may 
    (for values less than one) 
    make fitting more stable at the cost of efficiency.
    Can be set to "adaptive"; see 
    \code{MCMLE.adaptive.trustregion} and 
    \code{MCMLE.adaptive.epsilon}.
  }
  \item{MCMLE.adaptive.trustregion}{Maximum increase the algorithm
    will allow for the approximated loglikelihood at a given
    iteration when \code{MCMLE.steplength="adaptive"}.
  }
  \item{MCMLE.adaptive.epsilon}{convergence tolerance:  If
    the change in loglikelihood value 
    when \code{MCMLE.steplength="adaptive"} is smaller than this
    tolerance (and the adaptive steplength equals 1), 
    stop even if we haven't yet reached \code{MCMLE.maxit}
    iterations.
  }
  \item{MCMLE.sequential}{Logical:  If TRUE,
    the next iteration of the fit uses the last network 
    sampled as the starting network.  If FALSE, always use the 
    initially passed network.
    The results should be similar (stochastically), but the
    TRUE option may help if the \code{target.stats} 
    in the \code{\link{ergm}} function are far from the
    initial network.
  }
  \item{MCMLE.density.guard.min, MCMLE.density.guard}{
    A simple heuristic to stop optimization if it finds itself in an
    overly dense region, which usually indicates ERGM degeneracy: if the
    sampler encounters a network configuration that has more than
    \code{MCMLE.density.guard.min} edges and whose number of edges is
    exceeds the observed network by more than
    \code{MCMLE.density.guard}, the optimization process will be stopped
    with an error.
  }
  \item{SA.phase1_n}{Number of MCMC samples to draw
    in Phase 1 of the stochastic approximation algorithm.
    Defaults to 7 plus 3 times the number of terms in the model.
    See Snijders (2002) for details.
  }
  \item{SA.initial_gain}{Initial gain to
    Phase 2 of the stochastic approximation algorithm.
    See Snijders (2002) for details.
  }
  \item{SA.nsubphases}{Number of sub-phases 
    in Phase 2 of the stochastic approximation algorithm.
    Defaults to \code{MCMLE.maxit}.
    See Snijders (2002) for details.
  }
  \item{SA.niterations}{Number of MCMC samples to draw
    in Phase 2 of the stochastic approximation algorithm.
    Defaults to 7 plus the number of terms in the model.
    See Snijders (2002) for details.
  }
  \item{SA.phase3_n}{ample size for the MCMC sample
    in Phase 3 of the stochastic approximation algorithm.
    See Snijders (2002) for details.
  }
  \item{SA.trustregion}{Not yet documented.
  }

  \item{RM.phase1n_base}{The Robbins-Monro control
    parameters are not yet documented.
  }
  \item{RM.phase2n_base}{Not yet documented.
  }
  \item{RM.phase2sub}{Not yet documented.
  }
  \item{RM.init_gain}{Not yet documented.
  }
  \item{RM.phase3n}{Not yet documented.
  }
  \item{Step.MCMC.samplesize}{MCMC sample size for the preliminary steps of the
    "Stepping" method of optimization.  This is usually chosen to be smaller
    than the final MCMC sample size (which equals \code{MCMC.samplesize}).
    See Hummel et al. (2012) for details.
  }
  \item{Step.maxit}{Maximum number of iterations (steps) allowed by the
    "Stepping" method.  
  }
  \item{Step.gridsize}{Integer \eqn{N} such that the "Stepping" style of
    optimization chooses a step length equal to the largest possible multiple
    of \eqn{1/N}.
    See Hummel et al. (2012) for details.
  }
  \item{loglik.control}{See \code{\link{control.ergm.bridge}}
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
  \item{\dots}{Additional arguments, passed to other functions    
    This argument is helpful because it collects any control 
    parameters that have been deprecated; a warning message is
    printed in case of deprecated arguments.
  }
}
\value{
  A list with arguments as components.
}
\seealso{ \code{\link{ergm}}. The \code{\link{control.simulate}} 
  function performs a 
  similar function for
  \code{\link{simulate.ergm}}; 
  \code{\link{control.gof}} performs a
  similar function for \code{\link{gof}}.   
}
\references{
  \itemize{
    \item
    Snijders, T.A.B. (2002),
    Markov Chain Monte Carlo Estimation of Exponential Random Graph Models.
    Journal of Social Structure.
    Available from 
    \url{http://www.cmu.edu/joss/content/articles/volume3/Snijders.pdf}.
    
    \item Firth (1993),
    Bias Reduction in Maximum Likelihood Estimates.
    Biometrika,
    80: 27-38.

    \item Hunter, D. R. and M. S. Handcock (2006), Inference in curved
    exponential family models for networks. Journal of Computational
    and Graphical Statistics, 15: 565-583.
    
    \item Hummel, R. M., Hunter, D. R., and Handcock, M. S. (2012),
    Improving Simulation-Based Algorithms for Fitting ERGMs,
    Journal of Computational and Graphical Statistics, to appear.
  }
}
\keyword{ models }

