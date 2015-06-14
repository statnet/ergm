%  File man/control.ergm.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2015 Statnet Commons
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
                         SAN.burnin=MCMC.burnin * SAN.burnin.times,
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
             MCMLE.steplength.min=0.0001,

             SA.phase1_n=NULL, 
             SA.initial_gain=NULL, 
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
    \code{"zeros"} for a \code{0} vector of appropriate length.
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
  \item{MCMC.prop.weights, obs.MCMC.prop.weights}{Specifies the proposal distribution used in the MCMC
    Metropolis-Hastings algorithm.  Possible choices depending on
    selected \code{reference} and \code{constraints} arguments of the \code{\link{ergm}} function,
    but often include \code{"TNT"} and \code{"random"}, and the \code{"default"} is to use the one with the
    highest priority available.
    
    The \code{TNT} (tie / no tie) option puts roughly equal weight on selecting a
    dyad with or without a tie as a candidate for toggling, whereas the
    \code{random} option puts equal weight on all possible dyads, though the
    interpretation of \code{random} may change according to the constraints in
    place.  When no constraints are in place, the default is TNT, which 
    appears to improve Markov chain mixing particularly for networks with a
    low edge density, as is typical of many realistic social networks.

    \code{obs.MCMC.prop.weights}, if given separately, specifies the
    weights to be used for the constrained MCMC when missing dyads are
    present, defaulting to the same as \code{MCMC.prop.weights}.
  }
  \item{MCMC.prop.args, obs.MCMC.prop.args}{An alternative, direct way
    of specifying additional arguments to proposal. \code{obs.MCMC.prop.args}, if given separately, specifies the
    weights to be used for the constrained MCMC when missing dyads are
    present, defaulting to the same as \code{MCMC.prop.args}.
  }
  \item{MCMC.interval}{Number of proposals between sampled statistics. Increasing interval will reduces the autocorrelation in the sample, and may increase the precision in estimates by reducing MCMC error, at the expense of time. Set the interval higher for larger networks.
  }
  \item{MCMC.burnin}{Number of proposals before any MCMC sampling
    is done. It typically is set to a fairly large number.
  }
  \item{MCMC.samplesize}{Number of network statistics,
    randomly drawn from a given distribution on the set of all networks,
    returned by the Metropolis-Hastings algorithm. 
    Increasing sample size may increase the precision in the estimates by reducing MCMC error,
    at the expense of time. Set it higher for larger networks, or when using parallel functionality.
  }  
  \item{MCMLE.effectiveSize, MCMC.effectiveSize, MCMC.effectiveSize.damp, MCMC.effectiveSize.maxruns,
  MCMC.effectiveSize.base, MCMC.effectiveSize.points, MCMC.effectiveSize.order}{
    Set \code{MCMLE.effectiveSize} to non-NULL value to adaptively
             determine the burn-in and the MCMC length needed
    to get the specified effective size using the method of Sahlin (2011); 50 is a reasonable value.  This
             feature is in experimental status until we verify the
             coverage of the standard errors.

	     
  }  
  \item{MCMC.return.stats}{Logical:  If TRUE, return the matrix
    of MCMC-sampled network statistics.  This matrix should have
    \code{MCMC.samplesize} rows.
    This matrix can be 
    used directly by the \code{coda} package to assess MCMC
    convergence.
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
  \item{MCMC.addto.se}{Whether to add the standard errors induced by the
    MCMC algorithm to the estimates' standard errors.
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
  \item{SAN.burnin.times}{Multiplier for \code{SAN.burnin} relative to \code{MCMC.burnin}. This lets one control the amount of SAN burn-in (arguably, the most important of SAN parameters) without overriding the other SAN.control defaults. }
  \item{SAN.control}{Control arguments to \code{\link{san}}.
    See \code{\link{control.san}} for details.
  }

  \item{MCMLE.termination}{The criterion used for terminating MCMLE estimation: 
   \itemize{
      \item{"Hummel"}
      {Terminate when the Hummel step length
    is 1 for two consecutive iterations. For the last iteration, the sample size is
    boosted by a factor of \code{MCMLE.last.boost}. See Hummel et. al. (2012).

    Note that this criterion is incompatible with
   \code{MCMLE.steplength} \eqn{\ne} 1 or \code{MCMLE.steplength.margin}
   \eqn{=} \code{NULL}.}
      \item{"Hotelling"}
      {After every MCMC sample, an autocorrelation-adjusted Hotelling's T^2
    test for equality of MCMC-simulated network statistics to observed
    is conducted, and if its P-value exceeds \code{MCMLE.conv.min.pval}, the
    estimation is considered to have converged and finishes. This was the default option
    in \code{ergm} version 3.1.}
      \item{"precision"}
      {Terminate when the estimated loss in estimating precision due to
	using MCMC standard errors is below the precision bound specified by
      \code{MCMLE.MCMC.precision}, and the Hummel step length is 1 for two consecutive
      iterations. See \code{MCMLE.MCMC.precision} for details. This
   feature is in experimental status until we verify the coverage of the
   standard errors.

   Note that this criterion is incompatible with
   \eqn{\code{MCMLE.steplength}\ne 1} or \eqn{\code{MCMLE.steplength.margin}=\code{NULL}}.}
      \item{"none"}
      {Stop after \code{MCMLE.maxit} iterations.}
  }
  }
  \item{MCMLE.maxit}{Maximum number of times the parameter for the MCMC
    should be updated by maximizing the MCMC likelihood. At each
    step the parameter is changed to the values that maximizes
    the MCMC likelihood based on the current sample.
  }

  \item{MCMLE.conv.min.pval}{ The P-value used in the Hotelling test for early termination.
    }
  
  \item{MCMLE.NR.maxit, MCMLE.NR.reltol}{The method, maximum number of iterations and relative tolerance to use within the \code{optim} rountine
  in the MLE optimization. Note that by default, ergm uses \code{trust}, and falls back
  to \code{optim} only when \code{trust} fails.
  }  
  
  \item{obs.MCMC.samplesize, obs.MCMC.burnin, obs.MCMC.interval, obs.MCMC.burnin.min}{
    Sample size, burnin, and interval parameters for the MCMC sampling used
    when unobserved data are present in the estimation routine.    
  }
  \item{MCMLE.check.degeneracy}{Logical:  If TRUE, employ
    a check for model degeneracy.
  }
  \item{MCMLE.MCMC.precision, MCMLE.MCMC.max.ESS.frac}{\code{MCMLE.MCMC.precision} is a vector of upper bounds on the standard errors induced by the MCMC algorithm, expressed as a percentage of the total standard error. The MCMLE algorithm will terminate when the MCMC standard errors are below the
    precision bound, and the Hummel step length is 1 for two consecutive iterations. This is an experimental feature.
    
    If effective sample size is used (see \code{MCMC.effectiveSize}), then ergm may increase the target ESS
    to reduce the MCMC standard error.
  }
  \item{MCMLE.metric}{Method to calculate the loglikelihood approximation.
    See Hummel et al (2010)
    for an explanation of "lognormal" and "naive".
  }
  \item{MCMLE.method}{Deprecated. By default, ergm uses \code{trust}, and falls back
  to \code{optim} with Nelder-Mead method when \code{trust} fails.}
  
  \item{MCMLE.trustregion}{Maximum increase the algorithm will
    allow for the approximated likelihood at a given iteration.
    See Snijders (2002) for details.
  }
  \item{MCMLE.dampening}{(logical) Should likelihood dampening be used?}
  \item{MCMLE.dampening.min.ess}{The effective sample size below which dampening is used.}
  \item{MCMLE.dampening.level}{The proportional distance from boundary of the convex hull move.}
  \item{MCMLE.steplength.margin}{
    The extra margin required for a Hummel step to count as being inside 
    the convex hull of the sample. 
    Set this to 0 if the step length gets stuck at the same value over
    several iteraions. Set it to \code{NULL} to use fixed step
    length. Note that this parameter is required to be non-\code{NULL}
    for MCMLE termination using Hummel or precision criteria.
  }
  \item{MCMLE.steplength}{Multiplier for step length, which may 
    (for values less than one) 
    make fitting more stable at the cost of computational efficiency.
    Can be set to "adaptive"; see 
    \code{MCMLE.adaptive.trustregion}.

    If \code{MCMLE.steplength.margin} is not \code{NULL}, the step
    length will be set using the algorithm of Hummel et al. (2010). In
    that case, it will serve as the maximum step length
    considered. However, setting it to anything other than 1 will
    preclude using Hummel or precision as termination criteria.
  }
  \item{MCMLE.adaptive.trustregion}{Maximum increase the algorithm
    will allow for the approximated loglikelihood at a given
    iteration when \code{MCMLE.steplength="adaptive"}.
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
  \item{MCMLE.last.boost}{For the Hummel termination criterion, increase
  the MCMC sample size of the last iteration by this factor.}
  \item{MCMLE.Hummel.esteq}{For curved ERGMs, should the estimating function values be used to compute the Hummel step length? This allows the Hummel stepping algorithm converge when some sufficient statistics are at 0.}
  \item{MCMLE.steplength.min}{Stops MCMLE estimation when the step length gets stuck below this minimum value.}
  
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
  \item{SA.phase3_n}{Sample size for the MCMC sample
    in Phase 3 of the stochastic approximation algorithm.
    See Snijders (2002) for details.
  }
  \item{SA.trustregion}{The trust region parameter for the likelihood functions, used in the stochastic approximation algorithm.
  }

  \item{RM.phase1n_base, RM.phase2n_base, RM.phase2sub, RM.init_gain, RM.phase3n}{The Robbins-Monro control
    parameters are not yet documented.
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

  \item{CD.nsteps, CD.multiplicity}{
    Main settings for contrastive divergence to obtain initial values for the
    estimation: respectively, the number of Metropolis--Hastings steps to take before
    reverting to the starting value and the number of tentative
    proposals per step. Computational experiments indicate that
    increasing \code{CD.multiplicity} improves the estimate faster than
    increasing \code{CD.nsteps} --- up to a point --- but it also samples from the wrong
    distribution, in the sense that while as
    \code{CD.nsteps}\eqn{\rightarrow\infty}, the CD estimate approaches
    the MLE, this is not the case for \code{CD.multiplicity}.

    In practice, MPLE, when available, usually outperforms CD for even a very high
    \code{CD.nsteps} (which is, in turn, not very stable), so CD is useful primarily when
    MPLE is not available. This feature is to be considered
    experimental and in flux.

    The default values have been set experimentally, providing a
    reasonably stable, if not great, starting values.
  }
  \item{CD.nsteps.obs, CD.multiplicity.obs}{
    When there are missing dyads, \code{CD.nsteps} and
    \code{CD.multiplicity} must be set to a relatively high value, as
    the network passed is not necessarily a good start for
    CD. Therefore, these settings are in effect if there are missing dyads in the observed
    network, using a higher default number of steps.
  }
  \item{CD.maxit, CD.conv.min.pval, CD.NR.maxit, CD.NR.reltol,
    CD.metric, CD.method, CD.trustregion, CD.dampening, CD.dampening.min.ess,
    CD.dampening.level, CD.steplength.margin, CD.steplength, CD.adaptive.trustregion,
    CD.adaptive.epsilon}{
    Miscellaneous tuning parameters of the CD sampler and
    optimizer. These have the same meaning as their \code{MCMC.*}
    counterparts.

    Note that only the Hotelling's stopping criterion is implemented for
    CD.
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
    \code{"PSOCK"}. Defaults to using the \code{parallel} package with PSOCK clusters. See
    \code{\link{ergm-parallel}}
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
    Journal of Computational and Graphical Statistics, 21: 920-939.

    \item Kristoffer Sahlin. Estimating convergence of Markov chain
    Monte Carlo simulations. Master's Thesis. Stockholm University, 2011. \url{http://www2.math.su.se/matstat/reports/master/2011/rep2/report.pdf}

  }
}
\keyword{ models }

