#  File R/control.ergm.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################
#' Auxiliary function for fine-tuning ERGM fitting.
#' 
#' This function is only used within a call to the [ergm()] function.
#' See the Usage section in [ergm()] for details. Also see the
#' Details section about some of the interactions between its
#' arguments.
#'
#' Different estimation methods or components of estimation have
#' different efficient tuning parameters; and we generally want to use
#' the estimation controls to inform the simulation controls in
#' [control.simulate.ergm()]. To accomplish this, `control.ergm()` uses
#' method-specific controls, with the method identified by the prefix:
#' \describe{
#'
#' \item{`CD`}{Contrastive Divergence estimation \insertCite{Kr17u}{ergm}}
#'
#' \item{`MPLE`}{Maximum Pseudo-Likelihood Estimation \insertCite{StIk90p}{ergm}}
#'
#' \item{`MCMLE`}{Monte-Carlo MLE \insertCite{HuHa06i,HuHu12i}{ergm}}
#'
#' \item{`SA`}{Stochastic Approximation via Robbins--Monro \insertCite{RoMo51s,Sn02m}{ergm}}
#'
#' \item{`SAN`}{Simulated Annealing used when `target.stats` are specified for [ergm()]}
#'
#' \item{`obs`}{Missing data MLE \insertCite{HaGi10m}{ergm}}
#'
#' \item{`init`}{Affecting how initial parameter guesses are obtained}
#'
#' \item{`parallel`}{Affecting parallel processing}
#'
#' \item{`MCMC`}{Low-level MCMC simulation controls}
#'
#' }
#'
#' Corresponding `MCMC` controls will usually be overwritten by the
#' method-specific ones. After the estimation finishes, they will
#' contain the last MCMC parameters used.
#'
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
#' } Passing \code{control.ergm(init=coef(prev.fit))} can be used to
#' ``resume'' an uncoverged [ergm()] run, though `checkpoint` and
#' `resume` would be better under most circumstances.
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
#' Valid initial methods for a given reference are set by the
#' `InitErgmReference.*` function.
#' @param main.method One of "MCMLE" (default) or
#'   "Stochastic-Approximation".  Chooses the estimation method used
#'   to find the MLE.  \code{MCMLE} attempts to maximize an
#'   approximation to the log-likelihood function.
#'   \code{Stochastic-Approximation} are both stochastic approximation
#'   algorithms that try to solve the method of moments equation that
#'   yields the MLE in the case of an exponential family model. The
#'   direct use of the likelihood function has many theoretical
#'   advantages over stochastic approximation, but the choice will
#'   depend on the model and data being fit. See Handcock (2000) and
#'   Hunter and Handcock (2006) for details.
#'
#' @param force.main Logical: If TRUE, then force MCMC-based estimation method,
#' even if the exact MLE can be computed via maximum pseudolikelihood
#' estimation.
#' @param main.hessian Logical: If TRUE, then an approximate Hessian matrix is
#' used in the MCMC-based estimation method.
#'
#' @param MPLE.samplesize,init.MPLE.samplesize
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
#'   is used instead, All of these can be specified either as numbers or as
#'   `function(d,e)` taking the number of informative dyads and
#'   informative edges. Specifying or returning a larger number than
#'   the number of informative dyads is safe.
#'
#' @param MPLE.type One of `"glm"`, `"penalized"`, or
#' `"logitreg"`.  Chooses method of calculating MPLE.  `"glm"` is the
#' usual formal logistic regression called via [glm()], whereas
#' `"penalized"` uses the bias-reduced method of Firth (1993) as
#' originally implemented by Meinhard Ploner, Daniela Dunkler, Harry
#' Southworth, and Georg Heinze in the "logistf" package. `"logitreg"` is
#' an "in-house" implementation that is slower and probably less stable but
#' supports nonlinear logistic regression. It is invoked automatically when the
#' model has curved terms.
#' @param MPLE.maxit Maximum number of iterations for `"logitreg"`
#' implementation of MPLE.
#'
#' @param
#'   MPLE.nonident,MPLE.nonident.tol,MPLE.nonvar,MCMLE.nonident,MCMLE.nonident.tol,MCMLE.nonvar
#'   A rudimentary nonidentifiability/multicollinearity diagnostic. If
#'   `MPLE.nonident.tol > 0`, test the MPLE covariate matrix or the CD
#'   statistics matrix has linearly dependent columns via [QR
#'   decomposition][qr] with tolerance `MPLE.nonident.tol`. This is
#'   often (not always) indicative of a non-identifiable
#'   (multicollinear) model. If nonidentifiable, depending on
#'   `MPLE.nonident` issue a warning, an error, or a message
#'   specifying the potentially redundant statistics. Before the
#'   diagnostic is performed, covariates that do not vary (i.e.,
#'   all-zero columns) are dropped, with their handling controlled by
#'   `MPLE.nonvar`. The corresponding `MCMLE.*` arguments provide a
#'   similar diagnostic for the unconstrained MCMC sample's estimating
#'   functions.
#'
#' @param
#'   MPLE.covariance.method,MPLE.covariance.samplesize,MPLE.covariance.sim.burnin,MPLE.covariance.sim.interval
#'   Controls for estimating the MPLE covariance
#'   matrix. `MPLE.covariance method` determines the method, with
#'   `invHess` (the default) returning the covariance estimate
#'   obtained from the [glm()]. `Godambe` estimates the covariance
#'   matrix using the Godambe-matrix \insertCite{ScHu23c}{ergm}. This
#'   method is recommended for dyad-dependent models. Alternatively,
#'   `bootstrap` estimates standard deviations using a parametric
#'   bootstrapping approach \insertCite{@see @ScDe17e}{ergm}. The
#'   other parameters control, respectively, the number of networks to
#'   simulate, the MCMC burn-in, and the MCMC interval for `Godambe`
#'   and `bootstrap` methods.
#'
#' @param MPLE.check If `TRUE` (the default), perform the MPLE
#'   existence check described by \insertCite{ScHu23c;textual}{ergm}.
#'
#' @param MPLE.constraints.ignore If `TRUE`, MPLE will ignore all
#'   dyad-independent constraints except for those due to attributes
#'   missingness. This can be used to avert evaluating and storing the
#'   [`rlebdm`]s for very large networks except where absolutely
#'   necessary. Note that this can be very dangerous unless you know
#'   what you are doing.
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
#' @template control_MCMC_effectiveSize
#' 
#' @param
#'   MCMLE.effectiveSize,MCMLE.effectiveSize.interval_drop,MCMLE.burnin,MCMLE.interval,MCMLE.samplesize,MCMLE.samplesize.per_theta,MCMLE.samplesize.min
#'   Sets the corresponding `MCMC.*` parameters when
#'   `main.method="MCMLE"` (the default). Used because defaults may be
#'   different for different methods. `MCMLE.samplesize.per_theta`
#'   controls the MCMC sample size (not target effective size) as a
#'   function of the number of (curved) parameters in the model, and
#'   `MCMLE.samplesize.min` sets the minimum sample size regardless of
#'   their number.
#'
#' @param SA.burnin,SA.interval,SA.samplesize Sets the corresponding
#'   `MCMC.*` parameters when `main.method="Stochastic-Approximation"`.
#'
#' @param MCMC.return.stats Numeric: If positive, include an
#'   [`mcmc.list`] (two, if observational process was involved) of
#'   MCMC network statistics from the last iteration of network of the
#'   estimation. They will be thinned to have length of at most
#'   `MCMC.return.stats`. They are used for MCMC diagnostics.
#'
#' @param MCMC.runtime.traceplot Logical: If `TRUE`, plot traceplots of the MCMC
#' sample after every MCMC MLE iteration.
#' @template control_MCMC_maxedges
#' @param MCMC.addto.se Whether to add the standard errors induced by the MCMC
#' algorithm to the estimates' standard errors.
#' @param SAN.maxit When \code{target.stats} argument is passed to
#' [ergm()], the maximum number of attempts to use [san()]
#' to obtain a network with statistics close to those specified.
#' @param SAN.nsteps.times Multiplier for \code{SAN.nsteps} relative to
#' \code{MCMC.burnin}. This lets one control the amount of SAN burn-in
#' (arguably, the most important of SAN parameters) without overriding the
#' other `SAN` defaults.
#' @param SAN Control arguments to [san()].  See
#' [control.san()] for details.
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
#' 
#' * `"confidence"`: Performs an equivalence test to prove with level
#' of confidence \code{MCMLE.confidence} that the true value of the
#' deviation of the simulated mean value parameter from the observed
#' is within an ellipsoid defined by the inverse-variance-covariance
#' of the sufficient statistics multiplied by a scaling factor
#' `control$MCMLE.MCMC.precision` (which has a different default).
#' 
#' * `"none"` Stop after
#' \code{MCMLE.maxit} iterations.  
#' @param MCMLE.maxit Maximum number of times the parameter for the MCMC should
#' be updated by maximizing the MCMC likelihood. At each step the parameter is
#' changed to the values that maximizes the MCMC likelihood based on the
#' current sample.
#' @param MCMLE.conv.min.pval The P-value used in the Hotelling test for early
#' termination.
#' @param MCMLE.confidence The confidence level for declaring
#'   convergence for `"confidence"` methods.
#' @param MCMLE.min.depfac,MCMLE.sampsize.boost.pow When using adaptive MCMC effective size, and methods that increase the MCMC sample size, use `MCMLE.sampsize.boost.pow` as the power of the boost amount (relative to the boost of the target effective size), but ensure that sample size is no less than `MCMLE.min.depfac` times the target effective size.
#' @param MCMLE.confidence.boost The maximum increase factor in sample
#'   size (or target effective size, if enabled) when the
#'   `"confidence"` termination criterion is either not approaching
#'   the tolerance region or is unable to prove convergence.
#' @param MCMLE.confidence.boost.threshold,MCMLE.confidence.boost.lag Sample size or target effective size will be increaed if the distance from the tolerance region fails to decrease more than MCMLE.confidence.boost.threshold in this many successive iterations.
#' @param MCMLE.NR.maxit,MCMLE.NR.reltol The method, maximum number of
#' iterations and relative tolerance to use within the \code{optim} rountine in
#' the MLE optimization. Note that by default, ergm uses \code{trust}, and
#' falls back to \code{optim} only when \code{trust} fails.
#'
#' @param
#'   obs.MCMC.prop,obs.MCMC.prop.weights,obs.MCMC.prop.args,obs.MCMLE.effectiveSize,obs.MCMC.samplesize,obs.MCMC.burnin,obs.MCMC.interval,obs.MCMC.mul,obs.MCMC.samplesize.mul,obs.MCMC.burnin.mul,obs.MCMC.interval.mul,obs.MCMC.effectiveSize,obs.MCMLE.burnin,obs.MCMLE.interval,obs.MCMLE.samplesize,obs.MCMLE.samplesize.per_theta,obs.MCMLE.samplesize.min
#'   Corresponding MCMC parameters and settings used for the constrained sample when
#'   unobserved data are present in the estimation routine. By default, they are controlled by the `*.mul`
#'   parameters, as fractions of the corresponding settings for the
#'   unconstrained (standard) MCMC.
#'
#'   These can, in turn, be controlled by `obs.MCMC.mul`, which can be
#'   used to set the overal multiplier for the number of MCMC steps in
#'   the constrained sample; one half of its effect applies to the
#'   burn-in and interval and the other half to the total sample
#'   size. For example, for `obs.MCMC.mul=1/4` (the default),
#'   `obs.MCMC.samplesize` is set to \eqn{\sqrt{1/4}=1/2} that of
#'   `obs.MCMC.samplesize`, and `obs.MCMC.burnin` and
#'   `obs.MCMC.interval` are set to \eqn{\sqrt{1/4}=1/2} of their
#'   respective unconstrained sampling parameters. When
#'   `MCMC.effectiveSize` or `MCMLE.effectiveSize` are given, their
#'   corresponding `obs` parameters are set to them multiplied by
#'   `obs.MCMC.mul`.
#'
#'   Lastly, if `MCMLE.effectiveSize` is not NULL but
#'   `obs.MCMLE.effectiveSize` is, the constrained sample's target
#'   effective size is set adaptively to achieve a similar precision
#'   for the estimating functions as that achieved for the
#'   unconstrained.
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
#' @param MCMLE.steplength Multiplier for step length (on the
#'   mean-value parameter scale), which may (for values less than one)
#'   make fitting more stable at the cost of computational efficiency.
#'
#'   If \code{MCMLE.steplength.margin} is not \code{NULL}, the step
#'   length will be set using the algorithm of Hummel et
#'   al. (2010). In that case, it will serve as the maximum step
#'   length considered. However, setting it to anything other than 1
#'   will preclude using Hummel or precision as termination criteria.
#'
#' @param MCMLE.steplength.parallel Whether parallel multisection
#'   search (as opposed to a bisection search) for the Hummel step
#'   length should be used if running in multiple threads. Possible
#'   values (partially matched) are `"never"`, and
#'   (default) `"observational"` (i.e., when missing data MLE is
#'   used).
#'
#' @param MCMLE.steplength.solver The linear program solver to use for
#'   MCMLE step length calculation. Can be either `"glpk"` to use
#'   \CRANpkg{Rglpk} or `"lpsolve"` to use \CRANpkg{lpSolveAPI}.
#'   \CRANpkg{Rglpk} can be orders of magnitude faster, particularly
#'   for models with many parameters and with large sample sizes, so
#'   it is used where available; but it requires an external library
#'   to install under some operating systems, so fallback to
#'   \CRANpkg{lpSolveAPI} provided.
#'
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
#'
#' @param MCMLE.steplength.miss.sample In fitting the missing data
#'   MLE, the rules for step length become more complicated. In short,
#'   it is necessary for \emph{all} points in the constrained sample
#'   to be in the convex hull of the unconstrained (though they may be
#'   on the border); and it is necessary for their centroid to be in
#'   its interior. This requires checking a large number of points
#'   against whether they are in the convex hull, so to speed up the
#'   procedure, a sample is taken of the points most likely to be
#'   outside it.  This parameter specifies the sample size or a
#'   function of the unconstrained sample matrix to determine the
#'   sample size. If the parameter or the return value of the function
#'   has a length of 2, the first element is used as the sample size,
#'   and the second element is used in an early-termination heuristic,
#'   only continuing the tests until this many test points in a row
#'   did not yield a change in the step length.
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
#' @param SA.phase1_n A constant or a function of number of free
#'   parameters `q`, number of free canonical statistic `p`, and
#'   network size `n`, giving the number of MCMC samples to draw in
#'   Phase 1 of the stochastic approximation algorithm.  Defaults to
#'   \eqn{\max(200, 7+3p)}.  See Snijders (2002) for details.
#'
#' @param SA.initial_gain Initial gain to Phase 2 of the stochastic
#'   approximation algorithm. Defaults to 0.1. See Snijders (2002) for
#'   details.
#' @param SA.nsubphases Number of sub-phases in Phase 2 of the
#'   stochastic approximation algorithm.  Defaults to
#'   \code{MCMLE.maxit}.  See Snijders (2002) for details.
#'
#' @param SA.min_iterations,SA.max_iterations A constant or a function
#'   of number of free parameters `q`, number of free canonical
#'   statistic `p`, and network size `n`, giving the baseline numbers
#'   of iterations within each subphase of Phase 2 of the stochastic
#'   approximation algorithm. Default to \eqn{7+p} and \eqn{207+p},
#'   respectively.  See Snijders (2002) for details.
#'
#' @param SA.phase3_n Sample size for the MCMC sample in Phase 3 of
#'   the stochastic approximation algorithm.  See Snijders (2002) for
#'   details.
#'
#' @param CD.nsteps,CD.multiplicity Main settings for contrastive
#'   divergence to obtain initial values for the estimation:
#'   respectively, the number of Metropolis--Hastings steps to take
#'   before reverting to the starting value and the number of
#'   tentative proposals per step. Computational experiments indicate
#'   that increasing \code{CD.multiplicity} improves the estimate
#'   faster than increasing \code{CD.nsteps} --- up to a point --- but
#'   it also samples from the wrong distribution, in the sense that
#'   while as \code{CD.nsteps}\eqn{\rightarrow\infty}, the CD estimate
#'   approaches the MLE, this is not the case for
#'   \code{CD.multiplicity}.
#' 
#'   In practice, MPLE, when available, usually outperforms CD for
#'   even a very high \code{CD.nsteps} (which is, in turn, not very
#'   stable), so CD is useful primarily when MPLE is not
#'   available. This feature is to be considered experimental and in
#'   flux.
#' 
#'   The default values have been set experimentally, providing a
#'   reasonably stable, if not great, starting values.
#'
#' @param CD.nsteps.obs,CD.multiplicity.obs When there are missing dyads,
#' \code{CD.nsteps} and \code{CD.multiplicity} must be set to a relatively high
#' value, as the network passed is not necessarily a good start for CD.
#' Therefore, these settings are in effect if there are missing dyads in the
#' observed network, using a higher default number of steps.
#' 
#' @param CD.samplesize.per_theta,obs.CD.samplesize.per_theta,CD.maxit,CD.conv.min.pval,CD.NR.maxit,CD.NR.reltol,CD.metric,CD.method,CD.dampening,CD.dampening.min.ess,CD.dampening.level,CD.steplength.margin,CD.steplength,CD.steplength.parallel,CD.adaptive.epsilon,CD.steplength.esteq,CD.steplength.miss.sample,CD.steplength.min,CD.steplength.solver
#'   Miscellaneous tuning parameters of the CD sampler and
#'   optimizer. These have the same meaning as their `MCMLE.*` and
#'   `MCMC.*` counterparts.
#' 
#'   Note that only the Hotelling's stopping criterion is implemented
#'   for CD.
#' 
#' @param loglik See [control.ergm.bridge()]
#' @template term_options
#' @template control_MCMC_parallel
#' @template seed
#' @template control_MCMC_packagenames
#' @template control_dots
#'
#' @return A list with arguments as components.
#' @seealso [ergm()]. The [control.simulate()] function
#' performs a similar function for [simulate.ergm()];
#' [control.gof()] performs a similar function for [gof()].
#' @references \insertAllCited{}
#'
#' * Firth (1993), Bias Reduction in Maximum Likelihood Estimates.
#' Biometrika, 80: 27-38.
#' 
#' 
#' * Kristoffer Sahlin. Estimating convergence of Markov chain Monte Carlo
#' simulations. Master's Thesis. Stockholm University, 2011.
#' \url{https://www2.math.su.se/matstat/reports/master/2011/rep2/report.pdf}
#'
#' @keywords models
#' @export control.ergm
control.ergm<-function(drop=TRUE,

                       init=NULL,
                       init.method=NULL,
                       
                       main.method=c("MCMLE", "Stochastic-Approximation"),
                       force.main=FALSE,
                       main.hessian=TRUE,

                       checkpoint=NULL,
                       resume=NULL,

                       MPLE.samplesize=.Machine$integer.max,
                       init.MPLE.samplesize=function(d,e) max(sqrt(d),e,40)*8,
                       MPLE.type=c("glm", "penalized","logitreg"),
                       MPLE.maxit=10000,
                       MPLE.nonvar=c("warning","message","error"),
                       MPLE.nonident=c("warning","message","error"),
                       MPLE.nonident.tol=1e-10,
                       MPLE.covariance.samplesize =500,
                       MPLE.covariance.method ="invHess",
                       MPLE.covariance.sim.burnin = 1024,
                       MPLE.covariance.sim.interval = 1024,
                       MPLE.check = TRUE,
                       MPLE.constraints.ignore = FALSE,

                       MCMC.prop=trim_env(~sparse + .triadic),
                       MCMC.prop.weights="default", MCMC.prop.args=list(),
                       MCMC.interval=NULL,
                       MCMC.burnin=EVL(MCMC.interval*16),
                       MCMC.samplesize=NULL,
                       MCMC.effectiveSize=NULL,
                       MCMC.effectiveSize.damp=10,
                       MCMC.effectiveSize.maxruns=16,
                       MCMC.effectiveSize.burnin.pval=0.2,
                       MCMC.effectiveSize.burnin.min=0.05,
                       MCMC.effectiveSize.burnin.max=0.5,
                       MCMC.effectiveSize.burnin.nmin=16,
                       MCMC.effectiveSize.burnin.nmax=128,
                       MCMC.effectiveSize.burnin.PC=FALSE,
                       MCMC.effectiveSize.burnin.scl=32,
                       MCMC.effectiveSize.order.max=NULL,
                       MCMC.return.stats=2^12,
                       MCMC.runtime.traceplot=FALSE,
                       MCMC.maxedges=Inf,
                       MCMC.addto.se=TRUE,
                       MCMC.packagenames=c(),

                       SAN.maxit=4,
                       SAN.nsteps.times=8,
                       SAN=control.san(
                         term.options=term.options,
                         SAN.maxit=SAN.maxit,
                         SAN.prop=MCMC.prop,
                         SAN.prop.weights=MCMC.prop.weights,
                         SAN.prop.args=MCMC.prop.args,
                         
                         SAN.nsteps=EVL(MCMC.burnin,16384)*SAN.nsteps.times,
                         SAN.samplesize=EVL(MCMC.samplesize,1024),
                         SAN.packagenames=MCMC.packagenames,

                         parallel=parallel,
                         parallel.type=parallel.type,
                         parallel.version.check=parallel.version.check),
                       
                       MCMLE.termination=c("confidence", "Hummel", "Hotelling", "precision", "none"),
                       MCMLE.maxit=60,
                       MCMLE.conv.min.pval=0.5,
                       MCMLE.confidence=0.99,
                       MCMLE.confidence.boost=2,
                       MCMLE.confidence.boost.threshold=1,
                       MCMLE.confidence.boost.lag=4,
                       MCMLE.NR.maxit=100,
                       MCMLE.NR.reltol=sqrt(.Machine$double.eps),
                       obs.MCMC.mul=1/4,
                       obs.MCMC.samplesize.mul=sqrt(obs.MCMC.mul),
                       obs.MCMC.samplesize=EVL(round(MCMC.samplesize*obs.MCMC.samplesize.mul)),
                       obs.MCMC.effectiveSize=NVL3(MCMC.effectiveSize, .*obs.MCMC.mul),
                       obs.MCMC.interval.mul=sqrt(obs.MCMC.mul),
                       obs.MCMC.interval=EVL(round(MCMC.interval*obs.MCMC.interval.mul)),
                       obs.MCMC.burnin.mul=sqrt(obs.MCMC.mul),
                       obs.MCMC.burnin=EVL(round(MCMC.burnin*obs.MCMC.burnin.mul)),
                       obs.MCMC.prop=MCMC.prop, obs.MCMC.prop.weights=MCMC.prop.weights, obs.MCMC.prop.args=MCMC.prop.args,
                       obs.MCMC.impute.min_informative = function(nw) network.size(nw)/4,
                       obs.MCMC.impute.default_density = function(nw) 2/network.size(nw),

                       MCMLE.min.depfac=2,
                       MCMLE.sampsize.boost.pow=0.5,

                       MCMLE.MCMC.precision=if(startsWith("confidence", MCMLE.termination[1])) 0.1 else 0.005,
                       MCMLE.MCMC.max.ESS.frac=0.1,
                       MCMLE.metric=c("lognormal", "logtaylor",
                         "Median.Likelihood",
                         "EF.Likelihood", "naive"),
                       MCMLE.method=c("BFGS","Nelder-Mead"),
                       MCMLE.dampening=FALSE,
                       MCMLE.dampening.min.ess=20,
                       MCMLE.dampening.level=0.1,
                       MCMLE.steplength.margin=0.05,
                       MCMLE.steplength=NVL2(MCMLE.steplength.margin, 1, 0.5),
                       MCMLE.steplength.parallel=c("observational","never"),
                       MCMLE.sequential=TRUE,
                       MCMLE.density.guard.min=10000,
                       MCMLE.density.guard=exp(3),
                       MCMLE.effectiveSize=64,
                       obs.MCMLE.effectiveSize=NULL,
                       MCMLE.interval=1024,
                       MCMLE.burnin=MCMLE.interval*16,
                       MCMLE.samplesize.per_theta=32,
                       MCMLE.samplesize.min=256,
                       MCMLE.samplesize=NULL,
                       obs.MCMLE.samplesize.per_theta=round(MCMLE.samplesize.per_theta*obs.MCMC.samplesize.mul),
                       obs.MCMLE.samplesize.min=256,
                       obs.MCMLE.samplesize=NULL,
                       obs.MCMLE.interval=round(MCMLE.interval*obs.MCMC.interval.mul),
                       obs.MCMLE.burnin=round(MCMLE.burnin*obs.MCMC.burnin.mul),
                       MCMLE.steplength.solver=c("glpk","lpsolve"),
                       
                       MCMLE.last.boost=4,
                       MCMLE.steplength.esteq=TRUE, 
                       MCMLE.steplength.miss.sample=function(x1) c(max(ncol(rbind(x1))*2, 30), 10),
                       MCMLE.steplength.min=0.0001,
                       MCMLE.effectiveSize.interval_drop=2,
                       MCMLE.save_intermediates=NULL,
                       MCMLE.nonvar=c("message","warning","error"),
                       MCMLE.nonident=c("warning","message","error"),
                       MCMLE.nonident.tol=1e-10,

                       SA.phase1_n=function(q, ...) max(200, 7 + 3*q),
                       SA.initial_gain=0.1,
                       SA.nsubphases=4,
                       SA.min_iterations=function(q, ...) (7 + q),
                       SA.max_iterations=function(q, ...) (207 + q),
                       SA.phase3_n=1000,
                       SA.interval=1024,
                       SA.burnin=SA.interval*16,
                       SA.samplesize=1024,

                       CD.samplesize.per_theta=128,
                       obs.CD.samplesize.per_theta=128,
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
                       CD.dampening=FALSE,
                       CD.dampening.min.ess=20,
                       CD.dampening.level=0.1,
                       CD.steplength.margin=0.5,
                       CD.steplength=1,
                       CD.adaptive.epsilon=0.01,
                       CD.steplength.esteq=TRUE, 
                       CD.steplength.miss.sample=function(x1) ceiling(sqrt(ncol(rbind(x1)))),
                       CD.steplength.min=0.0001,
                       CD.steplength.parallel=c("observational","always","never"),
                       CD.steplength.solver=c("glpk","lpsolve"),
                       
                       loglik=control.logLik.ergm(),

                       term.options=NULL,

                       seed=NULL,
                       parallel=0,
                       parallel.type=NULL,
                       parallel.version.check=TRUE,
                       parallel.inherit.MT=FALSE,
                       
                       ...
                       ){
  old.controls <- list(SAN.control="SAN",
                       loglik.control="loglik",

                       CD.Hummel.esteq="CD.steplength.esteq",
                       CD.Hummel.miss.sample="CD.steplength.miss.sample",
                       MCMLE.Hummel.esteq="MCMLE.steplength.esteq",
                       MCMLE.Hummel.miss.sample="MCMLE.steplength.miss.sample",

                       mcmc.precision="MCMLE.MCMC.precision",
                       packagenames="MCMC.packagenames",
                       SAN.burnin.times="SAN.nsteps.times"
                       )

  for(trustarg in c("MCMLE.trustregion", "MCMLE.adaptive.trustregion",
                    "CD.trustregion", "CD.adaptive.trustregion",
                    "SA.trustregion"))
    old.controls[[trustarg]] <- list(action = warning, message = paste("The trust region mechanism has been obviated by step length", sQuote("*.steplen"), "and other mechanisms and has been removed."))
  old.controls[["MPLE.max.dyad.types"]] <- list(action = warning, message = paste("Argument", sQuote("MPLE.max.dyad.types"), " has been deprecated and will be removed in a future version."))

  match.arg.pars <- c("MPLE.type","MCMLE.metric","MCMLE.method","main.method",'MCMLE.termination',"CD.metric","CD.method","MCMLE.steplength.parallel","CD.steplength.parallel","MPLE.nonident","MPLE.nonvar","MCMLE.nonvar","MCMLE.nonident")

  control <- handle.controls("control.ergm", ...)

  if((control$MCMLE.steplength!=1 || is.null(control$MCMLE.steplength.margin)) && control$MCMLE.termination %in% c("Hummel", "precision"))
    stop("Hummel and precision-based termination require non-null MCMLE.steplength.margin and MCMLE.steplength = 1.")

  if(!is.null(control$checkpoint) && control$main.method!="MCMLE") stop("Only MCMLE supports checkpointing and resuming at this time.")

  set.control.class("control.ergm")
}


handle.control.toplevel<-function(myname, ...){
  myctrlname <- paste0("control.",myname)
  control.names <- ...names()[...names() %in% names(formals(get(myctrlname, mode="function")))]
  if(length(control.names)) stop("Argument(s) ", paste.and(sQuote(control.names)), " should be passed via control.",myname,"().")
}

SCALABLE_MCMC_CONTROLS <- c("MCMC.burnin", "MCMC.interval")
STATIC_MCMC_CONTROLS <- c("MCMC.samplesize", "MCMC.prop", "MCMC.prop.weights", "MCMC.prop.args", "MCMC.packagenames", "MCMC.maxedges", "term.options", "obs.MCMC.mul", "obs.MCMC.samplesize.mul", "obs.MCMC.samplesize", "obs.MCMC.interval.mul", "obs.MCMC.interval", "obs.MCMC.burnin.mul", "obs.MCMC.burnin", "obs.MCMC.prop", "obs.MCMC.prop.weights", "obs.MCMC.prop.args", "MCMC.batch")
ADAPTIVE_MCMC_CONTROLS <- c("MCMC.effectiveSize", "MCMC.effectiveSize.damp", "MCMC.effectiveSize.maxruns", "MCMC.effectiveSize.burnin.pval", "MCMC.effectiveSize.burnin.min", "MCMC.effectiveSize.burnin.max", "MCMC.effectiveSize.burnin.nmin", "MCMC.effectiveSize.burnin.nmax", "MCMC.effectiveSize.burnin.PC", "MCMC.effectiveSize.burnin.scl", "obs.MCMC.effectiveSize")
PARALLEL_MCMC_CONTROLS <- c("parallel","parallel.type","parallel.version.check")
OBS_MCMC_CONTROLS <- c("MCMC.base.samplesize", "MCMC.base.effectiveSize", "MCMC.samplesize", "MCMC.effectiveSize", "MCMC.interval", "MCMC.burnin")
MPLE_CONTROLS <- c("MPLE.samplesize", "MPLE.type", "MPLE.maxit", "drop")

remap_algorithm_MCMC_controls <- function(control, algorithm){
  CTRLS <- c(SCALABLE_MCMC_CONTROLS, STATIC_MCMC_CONTROLS, ADAPTIVE_MCMC_CONTROLS) %>% keep(startsWith,"MCMC.") %>% substr(6, 10000L)
  for(obs in c("", "obs.")){
    for(ctrl in CTRLS){
      dest <- paste0(obs, "MCMC.", ctrl)
      src <- paste0(obs, algorithm, ".", ctrl)
      if(length(control[[dest]])==0 && length(control[[src]])!=0) control[[dest]] <- control[[src]]
    }
  }
  control
}
