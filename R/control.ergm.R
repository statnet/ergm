#  File R/control.ergm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################
#' Auxiliary for Controlling ERGM Fitting
#' 
#' Auxiliary function as user interface for fine-tuning 'ergm' fitting.
#' 
#' This function is only used within a call to the [ergm()] function.
#' See the \code{usage} section in [ergm()] for details.
#' 
#' @templateVar MCMCType MCMC
#'
#' @param drop Logical: If TRUE, terms whose observed statistic values are at
#' the extremes of their possible ranges are dropped from the fit and their
#' corresponding parameter estimates are set to plus or minus infinity, as
#' appropriate.  This is done because maximum likelihood estimates cannot exist
#' when the vector of observed statistic lies on the boundary of the convex
#' hull of possible statistic values.
#' @param init numeric or \code{NA} vector equal in length to the number of
#' parameters in the model or \code{NULL} (the default); the initial values for
#' the estimation and coefficient offset terms. If \code{NULL} is passed, all
#' of the initial values are computed using the method specified by
#' \code{\link[=control.ergm]{control$init.method}}.  If a numeric vector is
#' given, the elements of the vector are interpreted as follows: \itemize{
#' \item Elements corresponding to terms enclosed in \code{offset()} are used as
#' the fixed offset coefficients. Note that offset coefficients alone can be
#' more conveniently specified using [ergm()] argument
#' \code{offset.coef}. If both \code{offset.coef} and \code{init} arguments are
#' given, values in \code{offset.coef} will take precedence.
#' 
#' \item Elements that do not correspond to offset terms and are not \code{NA}
#' are used as starting values in the estimation.
#' 
#' \item Initial values for the elements that are \code{NA} are fit using the
#' method specified by \code{\link[=control.ergm]{control$init.method}}.
#' 
#' } Passing \code{control.ergm(init=coef(prev.fit))} can be used to ``resume''
#' an uncoverged [ergm()] run, but see
#' \code{\link{enformulate.curved}}.
#' 
#' @param init.method A chatacter vector or \code{NULL}. The default
#'   method depends on the reference measure used. For the binary
#'   (\code{"Bernoulli"}) ERGMs, with dyad-independent constraints,
#'   it's maximum pseudo-likelihood estimation (MPLE). Other valid
#'   values include \code{"zeros"} for a \code{0} vector of
#'   appropriate length and \code{"CD"} for contrastive divergence. If
#'   passed explicitly, this setting overrides the reference's
#'   limitations.
#' 
#' Valid initial methods for a given reference are set by the `InitErgmReference.*` function.
#' @param main.method One of "MCMLE" (default),"Robbins-Monro",
#' "Stochastic-Approximation", or "Stepping".  Chooses the estimation method
#' used to find the MLE.  \code{MCMLE} attempts to maximize an approximation to
#' the log-likelihood function.  \code{Robbins-Monro} and
#' \code{Stochastic-Approximation} are both stochastic approximation algorithms
#' that try to solve the method of moments equation that yields the MLE in the
#' case of an exponential family model.  Another alternative is a partial
#' stepping algorithm (\code{Stepping}) as in Hummel et al. (2012).  The direct
#' use of the likelihood function has many theoretical advantages over
#' stochastic approximation, but the choice will depend on the model and data
#' being fit. See Handcock (2000) and Hunter and Handcock (2006) for details.
#' 
#' Note that in recent versions of ERGM, the enhancements of \code{Stepping}
#' have been folded into the default \code{MCMLE}, which is able to handle more
#' modeling scenarios.
#' @param force.main Logical: If TRUE, then force MCMC-based estimation method,
#' even if the exact MLE can be computed via maximum pseudolikelihood
#' estimation.
#' @param main.hessian Logical: If TRUE, then an approximate Hessian matrix is
#' used in the MCMC-based estimation method.
#'
#' @param MPLE.samplesize,init.MPLE.samplesize,MPLE.max.dyad.types
#'   These parameters control the maximum number of dyads (potential
#'   ties) that will be used by the MPLE to construct the predictor
#'   matrix for its logistic regression. In general, the algorithm
#'   visits dyads in a systematic sample that, if it does not hit one
#'   of these limits, will visit every informative dyad. If a limit is
#'   exceeded, case-control approximation to the likelihood,
#'   comprising all edges and those non-edges that have been visited
#'   by the algorithm before the limit was exceeded will be used.
#'
#'   `MPLE.samplesize` limits the number of dyads visited, unless the
#'   MPLE is being computed for the purpose of being the initial value
#'   for MCMC-based estimation, in which case `init.MPLE.samplesize`
#'   is used instead, `MPLE.max.dyad.types` limits the number of
#'   unique values of change statistic vectors that will be
#'   stored. All of these can be specified either as numbers or as
#'   `function(d,e)` taking the number of informative dyads and
#'   informative edges. Specifying or returning a larger number than
#'   the number of informative dyads is safe.
#'
#' @param MPLE.type One of `"glm"` or `"penalized"`.  Chooses method of calculating
#' MPLE.  `"glm"` is the usual formal logistic regression, whereas "penalized"
#' uses the bias-reduced method of Firth (1993) as originally implemented by
#' Meinhard Ploner, Daniela Dunkler, Harry Southworth, and Georg Heinze in the
#' "logistf" package.
#'
#' @param
#'   MPLE.nonident,MPLE.nonident.tol,MCMLE.nonident,MCMLE.nonident.tol
#'   A rudimentary nonidentifiability/multicollinearity diagnostic. If
#'   `MPLE.nonident.tol > 0`, test the MPLE covariate matrix or the CD
#'   statistics matrix has linearly dependent columns via [QR
#'   decomposition][qr] with tolerance `MPLE.nonident.tol`. This is
#'   often (not always) indicative of a non-identifiable
#'   (multicollinear) model. If nonidentifiable, depending on
#'   `MPLE.nonident` issue a warning, an error, or a message
#'   specifying the potentially redundant statistics. The
#'   corresponding `MCMLE.*` arguments provide a similar diagnostic
#'   for the unconstrained MCMC sample's estimating functions.
#'
#' @template control_MCMC_prop
#'
#' @param MCMC.interval Number of proposals between sampled statistics.
#' Increasing interval will reduces the autocorrelation in the sample, and may
#' increase the precision in estimates by reducing MCMC error, at the expense
#' of time. Set the interval higher for larger networks.
#' @param MCMC.burnin Number of proposals before any MCMC sampling is done. It
#' typically is set to a fairly large number.
#' @param MCMC.samplesize Number of network statistics, randomly drawn from a
#' given distribution on the set of all networks, returned by the
#' Metropolis-Hastings algorithm.  Increasing sample size may increase the
#' precision in the estimates by reducing MCMC error, at the expense of time.
#' Set it higher for larger networks, or when using parallel functionality.
#' @param
#' MCMLE.effectiveSize,MCMLE.effectiveSize.interval_drop,MCMC.effectiveSize,MCMC.effectiveSize.damp,MCMC.effectiveSize.maxruns,MCMC.effectiveSize.base,MCMC.effectiveSize.points,MCMC.effectiveSize.order
#' Set \code{MCMLE.effectiveSize} to non-NULL value to adaptively determine the
#' burn-in and the MCMC length needed to get the specified effective size using
#' the method of Sahlin (2011); 50 is a reasonable value.  This feature is in
#' experimental status until we verify the coverage of the standard errors.
#' 
#' @param MCMC.return.stats Logical: If TRUE, return the matrix of MCMC-sampled
#' network statistics.  This matrix should have \code{MCMC.samplesize} rows.
#' This matrix can be used directly by the \code{coda} package to assess MCMC
#' convergence.
#' @param MCMC.runtime.traceplot Logical: If `TRUE`, plot traceplots of the MCMC
#' sample after every MCMC MLE iteration.
#' @param MCMC.init.maxedges,MCMC.max.maxedges These parameters
#'   control how much space is allocated for storing edgelists for
#'   return at the end of MCMC sampling. Allocating more than needed
#'   wastes memory, so `MCMC.init.maxedges` is the initial amount
#'   allocated, but it will be incremented by a factor of 10 if
#'   exceeded during the simulation, up to `MCMC.max.maxedges`, at
#'   which point the process will stop with an error.
#' @param MCMC.addto.se Whether to add the standard errors induced by the MCMC
#' algorithm to the estimates' standard errors.
#' @param MCMC.compress Logical: If TRUE, the matrix of sample statistics
#' returned is compressed to the set of unique statistics with a column of
#' frequencies post-pended.
#' @param SAN.maxit When \code{target.stats} argument is passed to
#' [ergm()], the maximum number of attempts to use \code{\link{san}}
#' to obtain a network with statistics close to those specified.
#' @param SAN.nsteps.times Multiplier for \code{SAN.nsteps} relative to
#' \code{MCMC.burnin}. This lets one control the amount of SAN burn-in
#' (arguably, the most important of SAN parameters) without overriding the
#' other SAN.control defaults.
#' @param SAN.control Control arguments to \code{\link{san}}.  See
#' \code{\link{control.san}} for details.
#' @param MCMLE.termination The criterion used for terminating MCMLE
#' estimation:  
#' * `"Hummel"` Terminate when the Hummel step length is
#' 1 for two consecutive iterations. For the last iteration, the sample size is
#' boosted by a factor of \code{MCMLE.last.boost}. See Hummel et. al. (2012).
#' 
#' Note that this criterion is incompatible with \code{MCMLE.steplength}
#' \eqn{\ne} 1 or \code{MCMLE.steplength.margin} \eqn{=} \code{NULL}.
#' 
#' * `"Hotelling"` After every MCMC sample, an autocorrelation-adjusted
#' Hotelling's T^2 test for equality of MCMC-simulated network statistics to
#' observed is conducted, and if its P-value exceeds
#' \code{MCMLE.conv.min.pval}, the estimation is considered to have converged
#' and finishes. This was the default option in \code{ergm} version 3.1.
#' 
#' * `"precision"` Terminate when the estimated loss in estimating precision
#' due to using MCMC standard errors is below the precision bound specified by
#' \code{MCMLE.MCMC.precision}, and the Hummel step length is 1 for two
#' consecutive iterations. See \code{MCMLE.MCMC.precision} for details. This
#' feature is in experimental status until we verify the coverage of the
#' standard errors.
#' 
#' Note that this criterion is incompatible with
#' \eqn{\code{MCMLE.steplength}\ne 1} or
#' \eqn{\code{MCMLE.steplength.margin}=\code{NULL}}.  
#' * `"none"` Stop after
#' \code{MCMLE.maxit} iterations.  
#' @param MCMLE.maxit Maximum number of times the parameter for the MCMC should
#' be updated by maximizing the MCMC likelihood. At each step the parameter is
#' changed to the values that maximizes the MCMC likelihood based on the
#' current sample.
#' @param MCMLE.conv.min.pval The P-value used in the Hotelling test for early
#' termination.
#' @param MCMLE.NR.maxit,MCMLE.NR.reltol The method, maximum number of
#' iterations and relative tolerance to use within the \code{optim} rountine in
#' the MLE optimization. Note that by default, ergm uses \code{trust}, and
#' falls back to \code{optim} only when \code{trust} fails.
#' @param
#' obs.MCMC.prop.weights,obs.MCMC.prop.args,obs.MCMC.samplesize,obs.MCMC.burnin,obs.MCMC.interval,obs.MCMC.burnin.min
#' Corresponding MCMC parameters and settings used for the constrained sample when
#' unobserved data are present in the estimation routine.
#'
#' @param
#'   obs.MCMC.impute.min_informative,obs.MCMC.impute.default_density
#'   Controls for imputation of missing dyads for initializing MCMC
#'   sampling. If numeric, `obs.MCMC.impute.min_informative` specifies
#'   the minimum number dyads that need to be non-missing before
#'   sample network density is used as the imputation density. It can
#'   also be specified as a function that returns this
#'   value. `obs.MCMC.impute.default_density` similarly controls the
#'   imputation density when number of non-missing dyads is too low.
#' 
#' @param MCMLE.MCMC.precision,MCMLE.MCMC.max.ESS.frac
#' \code{MCMLE.MCMC.precision} is a vector of upper bounds on the standard
#' errors induced by the MCMC algorithm, expressed as a percentage of the total
#' standard error. The MCMLE algorithm will terminate when the MCMC standard
#' errors are below the precision bound, and the Hummel step length is 1 for
#' two consecutive iterations. This is an experimental feature.
#' 
#' If effective sample size is used (see \code{MCMC.effectiveSize}), then ergm
#' may increase the target ESS to reduce the MCMC standard error.
#' @param MCMLE.metric Method to calculate the loglikelihood approximation.
#' See Hummel et al (2010) for an explanation of "lognormal" and "naive".
#' @param MCMLE.method Deprecated. By default, ergm uses \code{trust}, and
#' falls back to \code{optim} with Nelder-Mead method when \code{trust} fails.
#' @param MCMLE.trustregion Maximum increase the algorithm will allow for the
#' approximated likelihood at a given iteration.  See Snijders (2002) for
#' details.
#' 
#' Note that not all metrics abide by it.
#' @param MCMLE.dampening (logical) Should likelihood dampening be used?
#' @param MCMLE.dampening.min.ess The effective sample size below which
#' dampening is used.
#' @param MCMLE.dampening.level The proportional distance from boundary of the
#' convex hull move.
#' @param MCMLE.steplength.margin The extra margin required for a Hummel step
#' to count as being inside the convex hull of the sample.  Set this to 0 if
#' the step length gets stuck at the same value over several iteraions. Set it
#' to \code{NULL} to use fixed step length. Note that this parameter is
#' required to be non-\code{NULL} for MCMLE termination using Hummel or
#' precision criteria.
#' @param MCMLE.steplength Multiplier for step length, which may (for values
#' less than one) make fitting more stable at the cost of computational
#' efficiency.  Can be set to "adaptive"; see
#' \code{MCMLE.adaptive.trustregion}.
#' 
#' If \code{MCMLE.steplength.margin} is not \code{NULL}, the step length will
#' be set using the algorithm of Hummel et al. (2010). In that case, it will
#' serve as the maximum step length considered. However, setting it to anything
#' other than 1 will preclude using Hummel or precision as termination
#' criteria.
#'
#' @param MCMLE.steplength.parallel Whether parallel multisection
#'   search (as opposed to a bisection search) for the Hummel step
#'   length should be used if running in multiple threads. Possible
#'   values (partially matched) are `"always"`, `"never"`, and
#'   (default) `"observational"` (i.e., when missing data MLE is
#'   used).
#'
#' @param MCMLE.adaptive.trustregion Maximum increase the algorithm will allow
#' for the approximated loglikelihood at a given iteration when
#' \code{MCMLE.steplength="adaptive"}.
#' @param MCMLE.sequential Logical: If TRUE, the next iteration of the fit uses
#' the last network sampled as the starting network.  If FALSE, always use the
#' initially passed network.  The results should be similar (stochastically),
#' but the TRUE option may help if the \code{target.stats} in the
#' [ergm()] function are far from the initial network.
#' @param MCMLE.density.guard.min,MCMLE.density.guard A simple heuristic to
#' stop optimization if it finds itself in an overly dense region, which
#' usually indicates ERGM degeneracy: if the sampler encounters a network
#' configuration that has more than \code{MCMLE.density.guard.min} edges and
#' whose number of edges is exceeds the observed network by more than
#' \code{MCMLE.density.guard}, the optimization process will be stopped with an
#' error.
#' @param MCMLE.last.boost For the Hummel termination criterion, increase the
#' MCMC sample size of the last iteration by this factor.
#' @param MCMLE.steplength.esteq For curved ERGMs, should the estimating function
#' values be used to compute the Hummel step length? This allows the Hummel
#' stepping algorithm converge when some sufficient statistics are at 0.
#' @param MCMLE.steplength.min Stops MCMLE estimation when the step length gets
#' stuck below this minimum value.
#' @param MCMLE.steplength.miss.sample In fitting the missing data MLE, the rules
#' for step length become more complicated. In short, it is necessary for
#' \emph{all} points in the constrained sample to be in the convex hull of the
#' unconstrained (though they may be on the border); and it is necessary for
#' their centroid to be in its interior. This requires checking a large number
#' of points against whether they are in the convex hull, so to speed up the
#' procedure, a sample is taken of the points most likely to be outside it.
#' This parameter specifies the sample size.
#' @param MCMLE.steplength.maxit Maximum number of iterations in searching for the
#' best step length.
#' 
#' @param checkpoint At the start of every iteration, save the state
#'   of the optimizer in a way that will allow it to be resumed. The
#'   name is passed through [sprintf()] with iteration number as the
#'   second argument. (For example, `checkpoint="step_%03d.RData"`
#'   will save to `step_001.RData`, `step_002.RData`, etc.)
#'
#' @param resume If given a file name of an `RData` file produced by
#'   `checkpoint`, the optimizer will attempt to resume after
#'   restoring the state. Control parameters from the saved state will
#'   be reused, except for those whose value passed via
#'   `control.ergm()` had change from the saved run. Note that if the
#'   network, the model, or some critical settings differ between
#'   runs, the results may be undefined.
#'
#' @param MCMLE.save_intermediates Every iteration, after MCMC
#'   sampling, save the MCMC sample and some miscellaneous information
#'   to a file with this name. This is mainly useful for diagnostics
#'   and debugging. The name is passed through [sprintf()] with
#'   iteration number as the second argument. (For example,
#'   `MCMLE.save_intermediates="step_%03d.RData"` will save to
#'   `step_001.RData`, `step_002.RData`, etc.)
#'
#' @param SA.phase1_n Number of MCMC samples to draw in Phase 1 of the
#' stochastic approximation algorithm.  Defaults to 7 plus 3 times the number
#' of terms in the model.  See Snijders (2002) for details.
#' @param SA.initial_gain Initial gain to Phase 2 of the stochastic
#' approximation algorithm.  See Snijders (2002) for details.
#' @param SA.nsubphases Number of sub-phases in Phase 2 of the stochastic
#' approximation algorithm.  Defaults to \code{MCMLE.maxit}.  See Snijders
#' (2002) for details.
#' @param SA.niterations Number of MCMC samples to draw in Phase 2 of the
#' stochastic approximation algorithm.  Defaults to 7 plus the number of terms
#' in the model.  See Snijders (2002) for details.
#' @param SA.phase3_n Sample size for the MCMC sample in Phase 3 of the
#' stochastic approximation algorithm.  See Snijders (2002) for details.
#' @param SA.trustregion The trust region parameter for the likelihood
#' functions, used in the stochastic approximation algorithm.
#' @param RM.phase1n_base,RM.phase2n_base,RM.phase2sub,RM.init_gain,RM.phase3n
#' The Robbins-Monro control parameters are not yet documented.
#' @param Step.MCMC.samplesize MCMC sample size for the preliminary steps of
#' the "Stepping" method of optimization.  This is usually chosen to be smaller
#' than the final MCMC sample size (which equals \code{MCMC.samplesize}).  See
#' Hummel et al. (2012) for details.
#' @param Step.maxit Maximum number of iterations (steps) allowed by the
#' "Stepping" method.
#' @param Step.gridsize Integer \eqn{N} such that the "Stepping" style of
#' optimization chooses a step length equal to the largest possible multiple of
#' \eqn{1/N}.  See Hummel et al. (2012) for details.
#' @param CD.nsteps,CD.multiplicity Main settings for contrastive divergence to
#' obtain initial values for the estimation: respectively, the number of
#' Metropolis--Hastings steps to take before reverting to the starting value
#' and the number of tentative proposals per step. Computational experiments
#' indicate that increasing \code{CD.multiplicity} improves the estimate faster
#' than increasing \code{CD.nsteps} --- up to a point --- but it also samples
#' from the wrong distribution, in the sense that while as
#' \code{CD.nsteps}\eqn{\rightarrow\infty}, the CD estimate approaches the MLE,
#' this is not the case for \code{CD.multiplicity}.
#' 
#' In practice, MPLE, when available, usually outperforms CD for even a very
#' high \code{CD.nsteps} (which is, in turn, not very stable), so CD is useful
#' primarily when MPLE is not available. This feature is to be considered
#' experimental and in flux.
#' 
#' The default values have been set experimentally, providing a reasonably
#' stable, if not great, starting values.
#' @param CD.nsteps.obs,CD.multiplicity.obs When there are missing dyads,
#' \code{CD.nsteps} and \code{CD.multiplicity} must be set to a relatively high
#' value, as the network passed is not necessarily a good start for CD.
#' Therefore, these settings are in effect if there are missing dyads in the
#' observed network, using a higher default number of steps.
#' 
#' @param CD.maxit,CD.conv.min.pval,CD.NR.maxit,CD.NR.reltol,CD.metric,CD.method,CD.trustregion,CD.dampening,CD.dampening.min.ess,CD.dampening.level,CD.steplength.margin,CD.steplength,CD.steplength.parallel,CD.adaptive.trustregion,CD.adaptive.epsilon,CD.steplength.esteq,CD.steplength.miss.sample,CD.steplength.maxit,CD.steplength.min
#'   Miscellaneous tuning parameters of the CD sampler and
#'   optimizer. These have the same meaning as their `MCMLE.*` and
#'   `MCMC.*` counterparts.
#' 
#'   Note that only the Hotelling's stopping criterion is implemented
#'   for CD.
#' 
#' @param loglik.control See \code{\link{control.ergm.bridge}}
#' @template term_options
#' @template control_MCMC_parallel
#' @template seed
#' @template control_MCMC_packagenames
#' @param \dots Additional arguments, passed to other functions This argument
#' is helpful because it collects any control parameters that have been
#' deprecated; a warning message is printed in case of deprecated arguments.
#' @return A list with arguments as components.
#' @seealso [ergm()]. The \code{\link{control.simulate}} function
#' performs a similar function for \code{\link{simulate.ergm}};
#' \code{\link{control.gof}} performs a similar function for \code{\link{gof}}.
#' @references \itemize{ 
#' * Snijders, T.A.B. (2002), Markov Chain Monte
#' Carlo Estimation of Exponential Random Graph Models.  Journal of Social
#' Structure.  Available from
#' \url{https://www.cmu.edu/joss/content/articles/volume3/Snijders.pdf}.
#' 
#' 
#' * Firth (1993), Bias Reduction in Maximum Likelihood Estimates.
#' Biometrika, 80: 27-38.
#' 
#' 
#' * Hunter, D. R. and M. S. Handcock (2006), Inference in curved
#' exponential family models for networks. Journal of Computational and
#' Graphical Statistics, 15: 565-583.
#' 
#' 
#' * Hummel, R. M., Hunter, D. R., and Handcock, M. S. (2012), Improving
#' Simulation-Based Algorithms for Fitting ERGMs, Journal of Computational and
#' Graphical Statistics, 21: 920-939.
#' 
#' 
#' * Kristoffer Sahlin. Estimating convergence of Markov chain Monte Carlo
#' simulations. Master's Thesis. Stockholm University, 2011.
#' \url{https://www2.math.su.se/matstat/reports/master/2011/rep2/report.pdf}
#' 
#' }
#' @keywords models
#' @export control.ergm
control.ergm<-function(drop=TRUE,

                       init=NULL,
                       init.method=NULL,
                       
                       main.method=c("MCMLE","Robbins-Monro",
                               "Stochastic-Approximation","Stepping"),
                       force.main=FALSE,
                       main.hessian=TRUE,

                       checkpoint=NULL,
                       resume=NULL,

                       MPLE.max.dyad.types=1e+6,
                       MPLE.samplesize=.Machine$integer.max,
                       init.MPLE.samplesize=function(d,e) max(sqrt(d),e,40)*8,
                       MPLE.type=c("glm", "penalized"),
                       MPLE.nonident=c("warning","message","error"),
                       MPLE.nonident.tol=1e-10,

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

                       SAN.maxit=4,
                       SAN.nsteps.times=8,
                       SAN.control=control.san(
                         term.options=term.options,
                         SAN.maxit=SAN.maxit,
                         SAN.prop.weights=MCMC.prop.weights,
                         SAN.prop.args=MCMC.prop.args,
                         SAN.init.maxedges=MCMC.init.maxedges,
                         SAN.max.maxedges=MCMC.max.maxedges,
                         
                         SAN.nsteps=MCMC.burnin*SAN.nsteps.times,
                         SAN.samplesize=MCMC.samplesize,
                         SAN.packagenames=MCMC.packagenames,

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
                       obs.MCMC.impute.min_informative = function(nw) network.size(nw)/4,
                       obs.MCMC.impute.default_density = function(nw) 2/network.size(nw),

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
                       MCMLE.steplength=NVL2(MCMLE.steplength.margin, 1, 0.5),
                       MCMLE.steplength.parallel=c("observational","always","never"),
                       MCMLE.adaptive.trustregion=3,
                       MCMLE.sequential=TRUE,
                       MCMLE.density.guard.min=10000,
                       MCMLE.density.guard=exp(3),
                       MCMLE.effectiveSize=NULL,
                       MCMLE.last.boost=4,
                       MCMLE.steplength.esteq=TRUE, 
                       MCMLE.steplength.miss.sample=100,
                       MCMLE.steplength.maxit=25, 
                       MCMLE.steplength.min=0.0001,
                       MCMLE.effectiveSize.interval_drop=2,
                       MCMLE.save_intermediates=NULL,
                       MCMLE.nonident=c("warning","message","error"),
                       MCMLE.nonident.tol=1e-10,

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
                       CD.steplength.esteq=TRUE, 
                       CD.steplength.miss.sample=100,
                       CD.steplength.maxit=25, 
                       CD.steplength.min=0.0001,
                       CD.steplength.parallel=c("observational","always","never"),
                       
                       loglik.control=control.logLik.ergm(),

                       term.options=NULL,

                       seed=NULL,
                       parallel=0,
                       parallel.type=NULL,
                       parallel.version.check=TRUE,
                       
                       ...
                       ){
  old.controls <- list(CD.Hummel.esteq="CD.steplength.esteq",
                       CD.Hummel.miss.sample="CD.steplength.miss.sample",
                       CD.Hummel.maxit="CD.steplength.maxit",
                       MCMLE.Hummel.esteq="MCMLE.steplength.esteq",
                       MCMLE.Hummel.miss.sample="MCMLE.steplength.miss.sample",
                       MCMLE.Hummel.maxit="MCMLE.steplength.maxit",

                       nr.maxit="MCMLE.NR.maxit",
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
                       packagenames="MCMC.packagenames",
                       SAN.burnin.times="SAN.nsteps.times"
                       )

  match.arg.pars <- c("MPLE.type","MCMLE.metric","MCMLE.method","main.method",'MCMLE.termination',"CD.metric","CD.method","MCMLE.steplength.parallel","CD.steplength.parallel","MPLE.nonident","MCMLE.nonident")
  
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

  if(!is.null(control$checkpoint) && control$main.method!="MCMLE") stop("Only MCMLE supports checkpointing and resuming at this time.")

  set.control.class("control.ergm")
}

control.toplevel<-function(..., myname= as.character(ult(sys.calls(), 2)[[1]])){
  myctrlname <- paste0("control.",myname) 
  control.names <- names(list(...))[names(list(...)) %in% names(formals(get(myctrlname, mode="function")))]
  if(length(control.names)) stop("Argument(s) ", paste.and(sQuote(control.names)), " should be passed via control.",myname,"().")
}
