#  File R/ergm.R in package ergm, part of the Statnet suite of packages for
#  network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################

#' Exponential-Family Random Graph Models
#'
#' [ergm()] is used to fit exponential-family random graph
#' models (ERGMs), in which
#' the probability of a given network, \eqn{y}, on a set of nodes is 
#' \eqn{h(y) \exp\{\eta(\theta) \cdot
#' g(y)\}/c(\theta)}, where
#' \eqn{h(y)} is the reference measure (usually \eqn{h(y)=1}),
#' \eqn{g(y)} is a vector of network statistics for \eqn{y},
#' \eqn{\eta(\theta)} is a natural parameter vector of the same 
#' length (with \eqn{\eta(\theta)=\theta} for most terms), and \eqn{c(\theta)} is the
#' normalizing constant for the distribution.
#' [ergm()] can return a maximum pseudo-likelihood
#' estimate, an approximate maximum likelihood estimate based on a Monte
#' Carlo scheme, or an approximate contrastive divergence estimate based
#' on a similar scheme.
#' (For an overview of the package \CRANpkg{ergm}, see \insertCite{HuHa08e;textual}{ergm} and \insertCite{KrHu23e;textual}{ergm}.)
#' \insertNoCite{AdHa07n,BeMo08p,Bu08sna,Bu08net,GoHa08s,GoKi09b,Ha03a,Ha03deg,HaGi10m,HaHu08s,HuHa06i,KaKr16s,Kr12e,MoHa08s,Sn02m}{ergm}
#'
#' @param formula An \R [`formula`], of the form
#'   \code{y ~ <model terms>}, where \code{y} is a
#'   [`network`] object or a matrix that can be
#'   coerced to a [`network`] object.  For the
#'   details on the possible \code{<model terms>}, see
#'   [`ergmTerm`] and \insertCite{MoHa08s}{ergm}
#'   for binary ERGM terms and \insertCite{Kr12e}{ergm} for valued ERGM terms
#'   (terms for weighted edges).  To create a
#'   [`network`] object in \R, use the
#'   \code{network()} function, then add nodal attributes to it using
#'   the \code{\%v\%} operator if necessary. Enclosing a model term in
#'   \code{offset()} fixes its value to one specified in
#'   \code{offset.coef}.  (A second argument---a logical or numeric
#'   index vector---can be used to select *which* of the parameters
#'   within the term are offsets.)
#'
#' @template response
#' @template reference
#' @template constraints
#'
#' @param obs.constraints A one-sided formula specifying one or more
#'   constraints or other modification \emph{in addition} to those
#'   specified by \code{constraints}, following the same syntax as the
#'   `constraints` argument.
#'
#'   This allows the domain of the integral in the numerator of the
#'   partially obseved network face-value likelihoods of
#'   \insertCite{HaGi10m}{ergm} and \insertCite{KaKr16s}{ergm} to be specified explicitly.
#'
#'   The default is to constrain the integral to only integrate over
#'   the missing dyads (if present), after incorporating constraints
#'   provided through the [`ergmlhs`] API.
#'
#'   It is also possible to specify a proposal function directly by
#'   passing a string with the function's name of the `obs.MCMC.prop`
#'   argument to the relevant control function. In that case,
#'   arguments to the proposal should be specified through the
#'   \code{obs.prop.args} argument to the relevant control function.
#'
#' @param offset.coef A vector of coefficients for the offset
#'   terms. Note that [`NaN`] elements are treated specially. See
#'   Skipping below. \matchnames{coefficient}
#'
#' @param target.stats vector of "observed network statistics,"
#' if these statistics are for some reason different than the
#' actual statistics of the network on the left-hand side of
#' \code{formula}.
#' Equivalently, this vector is the mean-value parameter values for the
#' model.  If this is given, the algorithm finds the natural
#' parameter values corresponding to these mean-value parameters.
#' If \code{NULL}, the mean-value parameters used are the observed
#' statistics of the network in the formula. \matchnames{statistic}
#'
#' @param eval.loglik {Logical:  For dyad-dependent models, if TRUE, use bridge
#' sampling to evaluate the log-likelihoood associated with the
#' fit. Has no effect for dyad-independent models.
#' Since bridge sampling takes additional time, setting to FALSE may
#' speed performance if likelihood values (and likelihood-based
#' values like AIC and BIC) are not needed. Can be set globally via `option(ergm.eval.loglik=...)`, which is set to `TRUE` when the package is loaded. (See [`options?ergm`][ergm-options].)
#' }
#' @param estimate {If "MPLE," then the maximum pseudolikelihood estimator
#' is returned.  If "MLE" (the default), then an approximate maximum likelihood
#' estimator is returned.  For certain models, the MPLE and MLE are equivalent,
#' in which case this argument is ignored.  (To force MCMC-based approximate
#' likelihood calculation even when the MLE and MPLE are the same, see the
#' \code{force.main} argument of [control.ergm()]. If "CD" (\emph{EXPERIMENTAL}),
#' the Monte-Carlo contrastive divergence estimate is returned. )
#' }
#'
#' @templateVar mycontrol control.ergm
#' @template control
#' @template verbose
#'
#' @param \dots Additional
#' arguments, to be passed to lower-level functions.
#'
#' @template basis
#'
#' @param newnetwork One of `"one"` (the default), `"all"`, or
#'   `"none"` (or, equivalently, `FALSE`), specifying whether the
#'   network(s) from the last iteration of the MCMC sampling should be
#'   returned as a part of the fit as a elements `newnetwork` and
#'   `newnetworks`. (See their entries in section Value below for
#'   details.) Partial matching is supported.
#'
#' @template MCMC-skipping
#'
#' @return [ergm()] returns an object of [`ergm`] that is a list
#'   consisting of the following elements: \item{coef}{The Monte Carlo
#'   maximum likelihood estimate of \eqn{\theta}, the vector of
#'   coefficients for the model parameters.}  \item{sample}{The
#'   \eqn{n\times p} matrix of network statistics, where \eqn{n} is
#'   the sample size and \eqn{p} is the number of network statistics
#'   specified in the model, generated by the last iteration of the
#'   MCMC-based likelihood maximization routine. These statistics are
#'   centered with respect to the observed statistics or
#'   `target.stats`, unless missing data MLE is used.}
#'   \item{sample.obs}{As \code{sample}, but for the constrained
#'   sample.}  \item{iterations}{The number of Newton-Raphson
#'   iterations required before convergence.}  \item{MCMCtheta}{The
#'   value of \eqn{\theta} used to produce the Markov chain Monte
#'   Carlo sample.  As long as the Markov chain mixes sufficiently
#'   well, \code{sample} is roughly a random sample from the
#'   distribution of network statistics specified by the model with
#'   the parameter equal to \code{MCMCtheta}.  If
#'   \code{estimate="MPLE"} then \code{MCMCtheta} equals the MPLE.}
#'   \item{loglikelihood}{The approximate change in log-likelihood in
#'   the last iteration.  The value is only approximate because it is
#'   estimated based on the MCMC random sample.}  \item{gradient}{The
#'   value of the gradient vector of the approximated loglikelihood
#'   function, evaluated at the maximizer.  This vector should be very
#'   close to zero.}  \item{covar}{Approximate covariance matrix for
#'   the MLE, based on the inverse Hessian of the approximated
#'   loglikelihood evaluated at the maximizer.}
#'   \item{failure}{Logical: Did the MCMC estimation fail?}
#'   \item{network}{Network passed on the left-hand side of
#'   `formula`. If `target.stats` are passed, it is replaced by the
#'   network returned by [san()].}  \item{newnetworks}{If argument
#'   `newnetwork` is `"all"`, a list of the final networks at the end
#'   of the MCMC simulation, one for each thread.}
#'   \item{newnetwork}{If argument `newnetwork` is `"one"` or `"all"`,
#'   the first (possibly only) element of \code{newnetworks}.}
#'   \item{coef.init}{The initial value of \eqn{\theta}.}
#'   \item{est.cov}{The covariance matrix of the model statistics in
#'   the final MCMC sample.}  \item{coef.hist, steplen.hist,
#'   stats.hist, stats.obs.hist}{ For the MCMLE method, the history of
#'   coefficients, Hummel step lengths, and average model statistics
#'   for each iteration..  } \item{control}{The control list passed to
#'   the call.}  \item{etamap}{The set of functions mapping the true
#'   parameter theta to the canonical parameter eta (irrelevant except
#'   in a curved exponential family model)} \item{formula}{The
#'   original [`formula`] passed to [ergm()].}
#'   \item{target.stats}{The target.stats used during estimation
#'   (passed through from the Arguments)} \item{target.esteq}{Used for
#'   curved models to preserve the target mean values of the curved
#'   terms. It is identical to target.stats for non-curved models.}
#'   \item{constraints}{Constraints used during estimation (passed
#'   through from the Arguments)} \item{reference}{The reference
#'   measure used during estimation (passed through from the
#'   Arguments)} \item{estimate}{The estimation method used (passed
#'   through from the Arguments).}  \item{offset}{vector of logical
#'   telling which model parameters are to be set at a fixed value
#'   (i.e., not estimated).}
#' 
#' \item{drop}{If \code{\link[=control.ergm]{control$drop}=TRUE}, a numeric vector indicating which terms were dropped due to to extreme values of the
#' corresponding statistics on the observed network, and how:
#' \describe{
#' \item{\code{0}}{The term was not dropped.}
#' \item{\code{-1}}{The term was at its minimum and the coefficient was fixed at
#' \code{-Inf}.}
#' \item{\code{+1}}{The term was at its maximum and the coefficient was fixed at
#' \code{+Inf}.}
#' }
#' }
#' 
#' \item{estimable}{A logical vector indicating which terms could not be
#' estimated due to a \code{constraints} constraint fixing that term at a
#' constant value.
#' }
#'
#' \item{info}{A list with miscellaneous information that would typically be accessed by the user via methods; in general, it should not be accessed directly. Current elements include: \describe{
#'
#' \item{`terms_dind`}{Logical indicator of whether the model terms are all dyad-independent.}
#'
#' \item{`space_dind`}{Logical indicator of whether the sample space (constraints) are all dyad-independent.}
#'
#' \item{`n_info_dyads`}{Number of \dQuote{informative} dyads: those that are observed (not missing) *and* not constrained by sample space constraints; one of the measures of sample size.}
#'
#' \item{`obs`}{Logical indicator of whether an observational (missing data) process was involved in estimation.}
#'
#' \item{`valued`}{Logical indicator of whether the model is valued.}
#'
#' }}
#' 
#' \item{null.lik}{Log-likelihood of the null model. Valid only for
#' unconstrained models.}
#' \item{mle.lik}{The approximate log-likelihood for the MLE.
#' The value is only approximate because it is estimated based 
#' on the MCMC random sample.}
#' 
#' @section Notes on model specification:
#' Although each of the statistics in a given model is a summary
#' statistic for the entire network, it is rarely necessary to
#' calculate statistics for an entire network
#' in a proposed Metropolis-Hastings step.
#' Thus, for example, if the triangle term is included in the model,
#' a census of all triangles in the observed network is never
#' taken; instead, only the change in the number of triangles
#' is recorded for each edge toggle.
#' 
#' In the implementation of [ergm()], the model is
#' initialized in \R, then all the model information is passed to a C
#' program that generates the sample of network statistics using MCMC.
#' This sample is then returned to \R, which then uses one of several
#' algorithms, selected by `main.method=` [control.ergm()] parameter
#' to update the estimate.
#' 
#' The mechanism for proposing new networks for the MCMC sampling
#' scheme, which is a Metropolis-Hastings algorithm, depends on 
#' two things:  The \code{constraints}, which define the set of possible
#' networks that could be proposed in a particular Markov chain step,
#' and the weights placed on these possible steps by the 
#' proposal distribution.  The former may be controlled using the
#' \code{constraints} argument described above.  The latter may
#' be controlled using the \code{prop.weights} argument to the
#' [control.ergm()] function.
#' 
#' The package is designed so that the user could conceivably add additional 
#' proposal types. 
#' 
#' @references \insertAllCited{}
#' 
#' @seealso [`network`], [`%v%`], [`%n%`], [`ergmTerm`], [`ergmMPLE`],
#' [summary.ergm()]
#' 
#' @examples
#' \donttest{
#' #
#' # load the Florentine marriage data matrix
#' #
#' data(flo)
#' #
#' # attach the sociomatrix for the Florentine marriage data
#' # This is not yet a network object.
#' #
#' flo
#' #
#' # Create a network object out of the adjacency matrix
#' #
#' flomarriage <- network(flo,directed=FALSE)
#' flomarriage
#' #
#' # print out the sociomatrix for the Florentine marriage data
#' #
#' flomarriage[,]
#' #
#' # create a vector indicating the wealth of each family (in thousands of lira) 
#' # and add it as a covariate to the network object
#' #
#' flomarriage %v% "wealth" <- c(10,36,27,146,55,44,20,8,42,103,48,49,10,48,32,3)
#' flomarriage
#' #
#' # create a plot of the social network
#' #
#' plot(flomarriage)
#' #
#' # now make the vertex size proportional to their wealth
#' #
#' plot(flomarriage, vertex.cex=flomarriage %v% "wealth" / 20, main="Marriage Ties")
#' #
#' # Use 'data(package = "ergm")' to list the data sets in a
#' #
#' data(package="ergm")
#' #
#' # Load a network object of the Florentine data
#' #
#' data(florentine)
#' #
#' # Fit a model where the propensity to form ties between
#' # families depends on the absolute difference in wealth
#' #
#' gest <- ergm(flomarriage ~ edges + absdiff("wealth"))
#' summary(gest)
#' #
#' # add terms for the propensity to form 2-stars and triangles
#' # of families 
#' #
#' gest <- ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle)
#' summary(gest)
#' 
#' # import synthetic network that looks like a molecule
#' data(molecule)
#' # Add a attribute to it to mimic the atomic type
#' molecule %v% "atomic type" <- c(1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,3)
#' #
#' # create a plot of the social network
#' # colored by atomic type
#' #
#' plot(molecule, vertex.col="atomic type",vertex.cex=3)
#' 
#' # measure tendency to match within each atomic type
#' gest <- ergm(molecule ~ edges + kstar(2) + triangle + nodematch("atomic type"))
#' summary(gest)
#' 
#' # compare it to differential homophily by atomic type
#' gest <- ergm(molecule ~ edges + kstar(2) + triangle
#'                         + nodematch("atomic type",diff=TRUE))
#' summary(gest)
#' }
#' @keywords models
#' @aliases is.ergm ergm.object
#' @export
ergm <- function(formula, response=NULL,
                 reference=~Bernoulli,
                 constraints=~.,
                 obs.constraints=~.-observed,
                 offset.coef=NULL,
                 target.stats=NULL,
                 eval.loglik=getOption("ergm.eval.loglik"),
                 estimate=c("MLE", "MPLE", "CD"),
                 control=control.ergm(),
                 verbose=FALSE,
                 ...,
                 basis=ergm.getnetwork(formula),
                 newnetwork=c("one", "all", "none")) {
  check_dots_used(error = unused_dots_warning)
  check.control.class("ergm", "ergm")
  handle.control.toplevel("ergm", ...)

  ergm_call <- match.call(ergm)

  reference <- trim_env_const_formula(reference)
  constraints <- trim_env_const_formula(constraints)
  obs.constraints <- trim_env_const_formula(obs.constraints)
  
  control$estimate <- estimate <- match.arg(estimate)
  newnetwork <-
    if(isFALSE(newnetwork)) "none"
    else if(is.character(newnetwork)) match.arg(newnetwork)

  if(estimate=="CD"){
    control$init.method <- "skip"
    control$main.method <- "CD"
    eval.loglik <- FALSE
  }

  if(estimate=="MPLE"){
    control$init.method <- "skip"
    control$main.method <- "MPLE"
  }
  
  if(!is.null(control$seed))  set.seed(as.integer(control$seed))
  if (verbose) message("Evaluating network in model.")
  
  nw <- basis
  ergm_preprocess_response(nw,response)

  proposalclass <- "c"

  if(!is(constraints, "ergm_proposal")){
    # Handle the observation process and other "automatic" constraints.
    tmp <- .handle.auto.constraints(nw, constraints, obs.constraints, target.stats)
    nw <- tmp$nw
    conterms.obs <- tmp$conterms.obs
    conterms <- tmp$conterms
  }else if(!is(obs.constraints, "ergm_proposal")){
    # Handle the observation process and other "automatic" constraints.
    tmp <- .handle.auto.constraints(nw, trim_env(~.), obs.constraints, target.stats)
    nw <- tmp$nw
    conterms.obs <- tmp$conterms.obs
    conterms <- tmp$conterms
  }
  
  if(!is(constraints, "ergm_proposal")){
    if (verbose) message("Initializing unconstrained Metropolis-Hastings proposal: ", appendLF=FALSE)
    
  ## FIXME: a more general framework is needed?
  if(is.valued(nw) && reference==trim_env(~Bernoulli)){
    warn(paste0("The default Bernoulli reference distribution operates in the binary (",sQuote("response=NULL"),") mode only. Did you specify the ",sQuote("reference")," argument?"))
  }
    
    proposal <- ergm_proposal(conterms, hints=control$MCMC.prop, weights=control$MCMC.prop.weights, control$MCMC.prop.args, nw, class=proposalclass,reference=reference, term.options=control$term.options)
  }else proposal <- constraints
  
  if (verbose) message(sQuote(paste0(proposal$pkgname,":MH_",proposal$name)),".")

  if(!is(obs.constraints, "ergm_proposal")){
    if(!is.null(conterms.obs)){
      if (verbose) message("Initializing constrained Metropolis-Hastings proposal: ", appendLF=FALSE)
      proposal.obs <- ergm_proposal(conterms.obs, hints=control$obs.MCMC.prop, weights=control$obs.MCMC.prop.weights, control$obs.MCMC.prop.args, nw, class=proposalclass, reference=reference, term.options=control$term.options)
    }else proposal.obs <- NULL
  }else proposal.obs <- obs.constraints

  if (verbose && !is.null(proposal.obs)) message(sQuote(paste0(proposal.obs$pkgname,":MH_",proposal.obs$name)),".")

  if (verbose) message("Initializing model...")
  model <- ergm_model(formula, nw, extra.aux = NVL3(proposal$auxiliaries,list(proposal=.)), term.options=control$term.options)
  proposal$aux.slots <- model$slots.extra.aux$proposal
  if (verbose) message("Model initialized.")

  model.obs <- NULL
  if(identical(proposal$auxiliaries, proposal.obs$auxiliaries)){
    ## Reuse auxiliaries from the unconstrained proposal if identical.
    proposal.obs$aux.slots <- model$slots.extra.aux$proposal
  }else if(!is.null(proposal.obs$auxiliaries)){
    if (verbose) message("Constrained proposal requires different auxiliaries: reinitializing model...")
    model.obs <- ergm_model(formula, nw, extra.aux = NVL3(proposal.obs$auxiliaries,list(proposal=.)), term.options=control$term.options)
    proposal.obs$aux.slots <- model.obs$slots.extra.aux$proposal
    if (verbose) message("Model reinitialized.")
  }

  info <- list(
    terms_dind = is.dyad.independent(model),
    space_dind = is.dyad.independent(proposal$arguments$constraints, proposal.obs$arguments$constraints),
    n_info_dyads = if(!control$MPLE.constraints.ignore) sum(as.rlebdm(proposal$arguments$constraints, proposal.obs$arguments$constraints, which="informative")) else NA,
    obs = !is.null(proposal.obs),
    valued = is.valued(nw)
  )
  info$MPLE_is_MLE <- (proposal$reference$name=="Bernoulli"
      && info$terms_dind
      && !control$force.main
      && info$space_dind)

  ## Construct approximate response network if target.stats are given.
  if(!is.null(target.stats)){
    target.stats <- match_names(target.stats, param_names(model, canonical=TRUE, offset=FALSE))
    
    # no need to pass the offset term's init to SAN
    san.control <- control$SAN
    
    if(verbose) message("Constructing an approximate response network.")
    ## If target.stats are given, overwrite the given network and formula
    ## with SAN-ed network and formula.
    if(control$SAN.maxit > 0){
      TARGET_STATS<-san(formula, target.stats=target.stats,
                reference=reference,
                constraints=constraints,
                control=san.control,
                only.last=TRUE,
                output="ergm_state",
                verbose=verbose,
                basis=nw,
                offset.coef=NVL(offset.coef,control$init[model$etamap$offsettheta]))
      if(verbose) message("Finished SAN run.")
    }else{
      TARGET_STATS <- nw
    }
    
    # From this point on, target.stats has NAs corresponding to the
    # offset terms.
    #
    # TODO: Only have target.stats contain non-offset terms'
    # statistics, and have the rest of the code handle it
    # intelligently.
    target.stats <- .embed.target.stats(model, target.stats)

    nw <- TARGET_STATS <- as.network(TARGET_STATS)
    #' @importFrom statnet.common nonsimp_update.formula
    formula<-nonsimp_update.formula(formula,TARGET_STATS~., from.new="TARGET_STATS")
  } else {
    if (network.edgecount(nw) == 0) warning("Network is empty and no target stats are specified.")
  }

  # The following kludge knocks out MPLE if the sample space
  # constraints are not dyad-independent. For example, ~observed
  # constraint is dyad-independent, while ~edges is not.
  #
  # TODO: Create a flexible and general framework to manage methods
  # for obtaining initial values.
  init.candidates <- c(proposal$reference$init_methods, "skip")
  if("MPLE" %in% init.candidates && !info$space_dind){
    init.candidates <- init.candidates %[.]% (. != "MPLE")
    if(verbose) message("MPLE cannot be used for this constraint structure.")
  }

  control$init.method <- ERRVL2(match.arg(control$init.method, init.candidates), {
    message("Specified initial parameter method ", sQuote(control$init.method), " is not in the list of candidates. Use at your own risk.")
    control$init.method
  })
  if(verbose && control$init.method!="skip") message(paste0("Using initial method '",control$init.method,"'."))

  if(is.curved(model) ||
     !all(model$etamap$mintheta==-Inf) ||
     !all(model$etamap$maxtheta==+Inf)){ # Curved or constrained model: use ergm.logitreg() rather than glm().
      control$MPLE.type <- "logitreg"
  }
  
  # If some control$init is specified...
  control$init <- match_names(NVL(control$init, structure(numeric(0), names=character(0))),
                              param_names(model, canonical=FALSE), NA)
 
  if(!is.null(offset.coef)){
    offset.coef <- match_names(offset.coef, param_names(model, canonical=FALSE, offset=TRUE))
    control$init[model$etamap$offsettheta]<-offset.coef
  }

  if(verbose){
    message("Initial parameters provided by caller:", appendLF = FALSE)
    init <- discard(control$init, is.na)
    if (length(init)) {
      message("")
      message_print(init)
    } else message(" None.")
    message("  number of free parameters: ", nparam(model, canonical = FALSE, offset = FALSE))
    message("  number of fixed parameters: ", nparam(model, canonical = FALSE, offset = TRUE))
  }

  # Make sure any offset elements are given in control$init.
  if(any(tmp <- is.NA(control$init) & model$etamap$offsettheta)) stop("The model contains offset terms whose parameter values have not been specified:", paste.and(sQuote(param_names(model)[tmp])), ".", sep="")

  # Make sure any NAN elements are offsets.
  if(any(tmp <- is.nan(control$init) & !model$etamap$offsettheta)) stop("The model contains non-offset terms whose parameter values are NAN:", paste.and(sQuote(param_names(model)[tmp])), ".", sep="")

  # Check if any terms are constrained to a constant and issue a warning.
  constrcheck <- ergm.checkconstraints.model(model, proposal, control$init)
  model <- constrcheck$model; control$init <- constrcheck$init
  
  # Check if any terms are at their extremes and handle them depending on control$drop.
  extremecheck <- ergm.checkextreme.model(model=model, nw=nw, init=control$init, target.stats=target.stats, drop=control$drop)
  model <- extremecheck$model; control$init <- extremecheck$init
  
  # MPLE is not supported for valued ERGMs.
  if(estimate=="MPLE"){
    if(is.valued(nw)) stop("Maximum Pseudo-Likelihood (MPLE) estimation for valued ERGM terms is not implemented at this time. You may want to use Contrastive Divergence by passing estimate='CD' instead.")
    if(!is.dyad.independent(proposal$arguments$constraints,
                            proposal.obs$arguments$constraints))
      message("Maximum Pseudo-Likelihood (MPLE) estimation for ERGMs with dyad-dependent constraints is not implemented at this time. You may want to use Contrastive Divergence by passing estimate='CD' instead.")
  }
  
  ## Run the fit.
  fit <- ergm.fit(nw, target.stats, model, model.obs, proposal, proposal.obs, info, control, verbose, ...)

  ## Process MCMC sample results.
  if(control$MCMC.return.stats == 0) fit$sample <- fit$sample.obs <- NULL
  else{ # Thin the chains.
    fit$sample <- NVL3(fit$sample, {
      w <- window(., thin = thin(.) * (et <- max(ceiling(niter(.) / control$MCMC.return.stats), 1)))
      structure(w, extra_thin = et)
    })
    fit$sample.obs <- NVL3(fit$sample.obs, {
      w <- window(., thin = thin(.) * (et <- max(ceiling(niter(.) / control$MCMC.return.stats), 1)))
      structure(w, extra_thin = et)
    })
  }

  ## Clear the new networks if not needed.
  if(newnetwork == "none") fit$newnetwork <- NULL
  if(newnetwork != "all") fit$newnetworks <- NULL

  ## Add miscellaneous information to the fit.
  fit$call <- ergm_call
  fit$network <- nw
  fit$ergm_version <- packageVersion("ergm")
  fit$info <- info
  fit$MPLE_is_MLE <- info$MPLE_is_MLE
  fit$drop <- if(control$drop) extremecheck$extremeval.theta
  fit$offset <- model$etamap$offsettheta
  fit$estimable <- constrcheck$estimable
  fit$etamap <- model$etamap
  fit$formula <- formula
  fit$target.stats <- suppressWarnings(na.omit(target.stats))
  fit$target.esteq <- suppressWarnings(na.omit(if(!is.null(target.stats)) ergm.estfun(rbind(target.stats), coef(fit), model)))
  fit$reference<-reference
  fit$constraints <- constraints
  fit$obs.constraints <- obs.constraints
  fit$estimate <- estimate
  fit$estimate.desc <- switch(estimate,
                              MPLE = if(info$MPLE_is_MLE) "Maximum Likelihood"
                                     else "Maximum Pseudolikelihood",
                              CD = "Contrastive Divergence",
                              MLE = paste0(if(NVL(fit$MCMCflag,FALSE)) # If not, it's just MLE.
                                             switch(control$main.method,
                                                    MCMLE = "Monte Carlo ",
                                                    `Stochastic-Approximation`="Stochastic Approximation "),
                                           "Maximum Likelihood"))

  NVL(fit$control) <- control

  class(fit) <- "ergm"

  ## Call the likelihood estimation.
  if(eval.loglik){
    message("Evaluating log-likelihood at the estimate. ", appendLF=FALSE)
    fit<-logLik(fit, add=TRUE, control=control$loglik, verbose=verbose)
    message("")
  }

  if (NVL(fit$MCMCflag,FALSE)) {
    message(paste(strwrap("This model was fit using MCMC.  To examine model diagnostics and check for degeneracy, use the mcmc.diagnostics() function."),collapse="\n"))
  }
  
  fit
}

ergm.fit <- function(nw, target.stats, model, model.obs, proposal, proposal.obs, info, control, verbose, ...){
  ## Short-circuit the optimization if all terms are either offsets or dropped.
  if(all(model$etamap$offsettheta)){
    ## Note that this cannot be overridden with control$force.main.
    message("All terms are either offsets or extreme values. No optimization is performed.")
    return(list(coefficients=control$init,
                iterations=0,
                loglikelihood=NA,
                gradient=rep(NA,length=length(control$init)),
                mple.lik=NA_real_
                ))
  }

  if(verbose) message("Fitting initial model.")

  ## Impute missing dyads and construct a state object.
  ##
  ## Note: We do not need to update nw.stats, because if we are in a
  ## situation where we are imputing dyads, the optimization is in the
  ## observational mode, and since both the constrained and the
  ## unconstrained samplers start from the same place, the initial
  ## stat shifts will be 0. target.stats and missing dyads are
  ## mutually exclusive, so target.stats will be set equal to
  ## network stats, causing this to happen.
  s <- single.impute.dyads(nw, constraints=proposal$arguments$constraints, constraints.obs=proposal.obs$arguments$constraints, min_informative = control$obs.MCMC.impute.min_informative, default_density = control$obs.MCMC.impute.default_density, output="ergm_state", verbose=verbose)

  ## statshift is the difference between the target.stats (if
  ## specified) and the statistics of the networks in the LHS of the
  ## formula or produced by SAN. If target.stats is not speficied
  ## explicitly, they are computed from this network, so
  ## statshift==0. To make target.stats play nicely with offsets, we
  ## set stat shifts to 0 where target.stats is NA (due to offset).
  nw.stats <- summary(model, s)
  statshift <- (nw.stats - NVL(target.stats, nw.stats)) |> replace(is.na, 0)

  s <- update(s, model=model, proposal=proposal, stats=statshift)
  s.obs <- NVL3(proposal.obs, update(s, model=NVL(model.obs,model), proposal=.))

  # If all other criteria for MPLE=MLE are met, _and_ SAN network
  # matches target.stats exactly, we can get away with MPLE.
  tgtmet <- NVL3(
    target.stats,
    isTRUE(all.equal(na.omit(.), nw.stats[!is.na(.)], check.attr = FALSE)),
    TRUE
  )

  if (!tgtmet) message("Unable to match target stats. Using MCMLE estimation.")
  if (control$estimate == "MLE"
      && info$MPLE_is_MLE && tgtmet
      && !control$force.main) {
    control$init.method <- "skip"
    control$main.method <- "MPLE"
  }

  ergm.getCluster(control, max(verbose-1,0)) # Set up parallel processing if necessary.
  
  initialfit <- ergm.initialfit(init=control$init,
                                s=s, s.obs=s.obs,
                                control=control,
                                verbose=max(verbose-1,0),
                                ...)

  ## Extract and process the initial value for the next stage:
  init <- coef(initialfit) |> replace(is.NA, 0)
  names(init) <- param_names(model, FALSE)

  c(
    switch(control$main.method,
           "MPLE" = ergm.mple(s, s.obs,
                              init=init,
                              control=control,
                              verbose=verbose, ...),

           "CD" = ergm.CD.fixed(.constrain_init(s$model, ifelse(is.NA(init),0,init)),
                      s, s.obs, control, verbose,...),

           "Stochastic-Approximation" = ergm.stocapprox(init, s, s.obs,
                                                        control=control,
                                                        verbose=verbose,
                                                        ...),

           "MCMLE" = ergm.MCMLE(init, s, s.obs,
                              control=control,
                              verbose=verbose,
                              ...),
                    
           stop("Method ", control$main.method, " is not implemented.")),
    list(nw.stats=nw.stats)
  )
}
