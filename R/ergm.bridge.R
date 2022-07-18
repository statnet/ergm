#  File R/ergm.bridge.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################

#' Bridge sampling to evaluate ERGM log-likelihoods and log-likelihood ratios
#' 
#' \code{ergm.bridge.llr} uses bridge sampling with geometric spacing to
#' estimate the difference between the log-likelihoods of two parameter vectors
#' for an ERGM via repeated calls to \code{\link{simulate.formula.ergm}}.
#' 
#' 
#' 
#' @param object A model formula. See \code{\link{ergm}} for details.
#' @template response
#' @param constraints,obs.constraints One-sided formulas specifying
#'   one or more constraints on the support of the distribution of the
#'   networks being simulated and on the observation process
#'   respectively. See the documentation for similar arguments for
#'   \code{\link{ergm}} for more information.
#' @param reference {A one-sided formula specifying the reference
#'   measure (\eqn{h(y)}) to be used.  (Defaults to
#'   \code{~Bernoulli}.)}
#' @param target.stats {A vector of sufficient statistics to be used
#'   in place of those of the network in the formula.}
#' @param from,to The initial and final parameter vectors.
#' @param basis An optional \code{\link[network]{network}} object to
#'   start the Markov chain.  If omitted, the default is the
#'   left-hand-side of the \code{object}.
#' @template verbose
#' @param \dots Further arguments to \code{ergm.bridge.llr} and
#'   \code{\link{simulate.formula.ergm}}.
#' @param llronly Logical: If TRUE, only the estiamted log-ratio will
#'   be returned by `ergm.bridge.llr`.
#'
#' @templateVar mycontrol control.ergm.bridge
#' @template control
#' @template verbose
#'
#' @param coef A vector of coefficients for the configuration of
#'   interest.
#' @param llkonly Whether only the estiamted log-likelihood should be
#'   returned by the `ergm.bridge.0.llk` and
#'   `ergm.bridge.dindstart.llk`.  (Defaults to TRUE.)
#' @return If `llronly=TRUE` or `llkonly=TRUE`, these functions return
#'   the scalar log-likelihood-ratio or the log-likelihood.
#'   Otherwise, they return a list with the following components:
#'
#'   \item{llr}{The estimated log-ratio.}
#'
#'   \item{llr.vcov}{The estimated variance of the log-ratio due to
#'   MCMC approximation.}
#'
#'   \item{llrs}{A list of lists (1 per attempt) of the estimated
#'   log-ratios for each of the \code{bridge.nsteps} bridges.}
#'
#'   \item{llrs.vcov}{A list of lists (1 per attempt) of the estimated
#'   variances of the estimated log-ratios for each of the
#'   \code{bridge.nsteps} bridges.}
#'
#'   \item{paths}{A list of lists (1 per attempt) with two elements:
#'   `theta`, a numeric matrix with `bridge.nsteps` rows, with each
#'   row being the respective bridge's parameter configuration; and
#'   `weight`, a vector of length `bridge.nsteps` containing its
#'   weight.}
#'
#'   \item{Dtheta.Du}{The gradient vector of the parameter values with
#'   respect to position of the bridge.}
#' @seealso \code{\link{simulate.formula.ergm}}
#' @references Hunter, D. R. and Handcock, M. S. (2006)
#'   \emph{Inference in curved exponential family models for
#'   networks}, Journal of Computational and Graphical Statistics.
#' @keywords model
#' @export
ergm.bridge.llr<-function(object, response=NULL, reference=~Bernoulli, constraints=~., from, to, obs.constraints=~.-observed, target.stats=NULL, basis=ergm.getnetwork(object), verbose=FALSE, ..., llronly=FALSE, control=control.ergm.bridge()){
  check.control.class("ergm.bridge", "ergm.bridge.llr")
  handle.control.toplevel("ergm.bridge", ...)

  if(!is.null(control$seed)) {set.seed(as.integer(control$seed))}

  # Set up a cluster, if not already.
  ergm.getCluster(control, verbose)

  message("Setting up bridge sampling...")
  
  ergm_preprocess_response(basis, response)

  ## Generate a path of n points shifted by shift in [-1/2, +1/2] as a
  ## fraction of the width of a bridge.
  ##
  ## The weight of each bridge is the size of the Voronoi partition
  ## for that bridge.
  mkpath <- function(n, shift = 0, reverse = FALSE) {
    stopifnot(shift >= -1/2, shift <= 1/2)
    u0 <- seq(from = 0 + 1 / 2 / n, to = 1 - 1 / 2 / n, length.out = n)
    u <- u0 + shift / n
    if (reverse) u <- rev(u)
    list(
      theta = t(rbind(sapply(u, function(u) cbind(to * u + from * (1 - u))))),
      u = u
    )
  }

  uweights <- function(u) {
    o <- order(u)
    u <- u[o]
    halfgap <- c(u[1], diff(u)/2, 1 - ult(u))
    (head(halfgap, -1) + tail(halfgap, -1))[order(o)]
  }

  # A low-discrepancy sequence: Kronecker Recurrence using inverse
  # Golden Ratio (0.6180...) step size.
  KR <- function(i, shift = 0) (shift + 0.61803398874989479 * i) %% 1
  # Ensure that shifts are in [-1/2, +1/2] and the first shift is 0.
  pathshift <- function(i) KR(i, shift = -KR(1) + 1/2) - 1/2

  # Determine whether an observation process is in effect.
  obs <- has.obs.constraints(basis, constraints, obs.constraints, target.stats)

  ## Control list constructor
  gen_control <- function(obs, burnin = c("first", "between")){
    control$MCMC.burnin <-
      switch(match.arg(burnin),
             first = if(obs) control$obs.MCMC.burnin
                     else control$MCMC.burnin,
             between = if(obs) control$obs.MCMC.burnin.between
                       else control$MCMC.burnin.between)
    control$MCMC.interval <-
      if(obs) ceiling(control$obs.MCMC.interval / control$bridge.nsteps)
      else ceiling(control$MCMC.interval / control$bridge.nsteps)

    #' @importFrom utils head
    modifyList(do.call(control.simulate.formula, control[intersect(names(control), head(names(formals(control.simulate.formula)), -1))]),
               list(MCMC.samplesize = if(obs) control$obs.MCMC.samplesize else control$MCMC.samplesize))
  }

  ## Obtain simulation setting arguments in terms of ergm_state.
  if(verbose) message("Initializing model and proposals...")
  sim_settings <- simulate(object, coef=from, nsim=1, reference=reference, constraints=list(constraints, obs.constraints), observational=FALSE, output="ergm_state", verbose=max(verbose-1,0), basis = basis, control=gen_control(FALSE, "first"), ..., return.args = "ergm_state")
  if(verbose) message("Model and proposals initialized.")
  state <- list(sim_settings$object)

  if(obs){
    if(verbose) message("Initializing constrained model and proposals...")
    sim_settings.obs <- simulate(object, coef=from, nsim=1, reference=reference, constraints=list(constraints, obs.constraints), observational=TRUE, output="ergm_state", verbose=max(verbose-1,0), basis = basis, control=gen_control(TRUE, "first"), ..., return.args = "ergm_state")
    if(verbose) message("Constrained model and proposals initialized.")
    state.obs <- list(sim_settings.obs$object)
  }

  ## Miscellaneous settings
  Dtheta.Du <- (to-from)[!state[[1]]$model$etamap$offsettheta]

  ## Handle target statistics, if passed.
  if(!is.null(target.stats)){
    if(nparam(as.ergm_model(state[[1]]), canonical=TRUE, offset=FALSE)!=length(target.stats)){
      stop("Incorrect length of the target.stats vector: should be ", nparam(as.ergm_model(state[[1]]), canonical=TRUE, offset=FALSE), " but is ",length(target.stats),". Note that offset() terms should *not* get target statistics.")
    }
    target.stats <- .align.target.stats.offset(as.ergm_model(state[[1]]), target.stats)
    if(any(as.ergm_model(state[[1]])$etamap$offsetmap)) warning("Using target.stats for a model with offset terms may produce an inaccurate estimate of the log-likelihood and derived quantities (deviance, AIC, BIC, etc.), because some of the target stats must be imputed.")
  }else target.stats <- summary(state[[1]])


  ## Helper function to calculate Dtheta.Du %*% Deta.Dtheta %*% g(y)
  llrsamp <- function(samp, theta){
    if(is.mcmc.list(samp)) lapply.mcmc.list(ergm.estfun(samp, theta, state[[1]]$model$etamap), `%*%`, Dtheta.Du)
    else sum(ergm.estfun(samp, theta, state[[1]]$model$etamap) * Dtheta.Du)
  }


  message("Using ", control$bridge.nsteps, " bridges: ", appendLF=FALSE)

  llr.hist <- list()
  vcov.llr.hist <- list()
  path.hist <- list()

  repeat{
    attempt <- length(path.hist) + 1
    # Bridge in reverse order on even-numbered attempts, if bidirectional bridging used.
    path <- mkpath(control$bridge.nsteps, shift = pathshift(attempt), reverse = control$bridge.bidirectional && attempt %% 2 == 0)
    llrs <- numeric(control$bridge.nsteps)
    vcov.llrs <- numeric(control$bridge.nsteps)

    for(i in seq_len(control$bridge.nsteps)){
      theta <- path$theta[i, ]
      if(verbose==0) message(i," ",appendLF=FALSE)
      if(verbose>0) message("Running theta=[",paste(format(theta),collapse=","),"].")

      ## First burn-in has to be longer, but those thereafter should be shorter if the bridges are closer together.
      z <- ergm_MCMC_sample(state, theta = theta, verbose = max(verbose - 1, 0),
                            control = gen_control(FALSE, if(i == 1 && (attempt == 1 || !control$bridge.bidirectional)) "first" else "between"))
      state <- z$networks
      samp <- llrsamp(z$stats, theta)
      vcov.llrs[i] <- c(ERRVL(try(spectrum0.mvar(samp)/(niter(samp)*nchain(samp)), silent=TRUE), 0))
      llrs[i] <- mean(as.matrix(samp))

      if(obs){
        z <- ergm_MCMC_sample(state.obs, theta = theta, verbose = max(verbose - 1, 0),
                              control = gen_control(TRUE, if(i == 1 && (attempt == 1 || !control$bridge.bidirectional)) "first" else "between"))
        state.obs <- z$networks
        samp <- llrsamp(z$stats, theta)
        vcov.llrs[i] <- vcov.llrs[i] + c(ERRVL(try(spectrum0.mvar(samp)/(niter(samp)*nchain(samp)), silent=TRUE), 0))
        llrs[i] <- llrs[i] - mean(as.matrix(samp))
      }else llrs[i] <- llrs[i] - llrsamp(target.stats, theta)
    }
    message(".")

    if(verbose) message("Bridge sampling finished. Collating...")

    llr.hist[[attempt]] <- llrs
    vcov.llr.hist[[attempt]] <- vcov.llrs
    path.hist[[attempt]] <- path

    w <- uweights(unlist(lapply(path.hist, `[[`, "u")))
    llr <- sum(unlist(llr.hist)*w)
    vcov.llr <- sum(unlist(vcov.llr.hist)*w^2)

    if(is.null(control$bridge.target.se) || vcov.llr <= control$bridge.target.se^2) break
    else message("Estimated standard error (", format(sqrt(vcov.llr)), ") above target (", format(control$bridge.target.se), "). Drawing additional samples.")
  }

  if(llronly) structure(llr, vcov=vcov.llr)
  else list(llr = llr, vcov.llr = vcov.llr, from = from, to = to, llrs = llr.hist, vcov.llrs = vcov.llr.hist, paths = path.hist, Dtheta.Du = Dtheta.Du)
}

#' @rdname ergm.bridge.llr
#'
#' @description \code{ergm.bridge.0.llk} is a convenience wrapper that
#'   returns the log-likelihood of configuration \eqn{\theta}
#'   \emph{relative to the reference measure}. That is, the
#'   configuration with \eqn{\theta=0} is defined as having log-likelihood of
#'   0.
#'
#' @return \code{ergm.bridge.0.llk} result list also includes an `llk`
#'   element, with the log-likelihood itself (with the reference
#'   distribution assumed to have likelihood 0).
#' 
#' 
#' @export
ergm.bridge.0.llk<-function(object, response=NULL, reference=~Bernoulli, coef, ..., llkonly=TRUE, control=control.ergm.bridge(), basis=ergm.getnetwork(object)){
  check.control.class("ergm.bridge", "ergm.bridge.0.llk")
  handle.control.toplevel("ergm.bridge", ...)
  ergm_preprocess_response(basis, response)
  br<-ergm.bridge.llr(object, from=rep(0,length(coef)), to=coef, reference=reference, control=control, ..., basis=basis)
  if(llkonly) br$llr
  else c(br,llk=br$llr)
}

#' @rdname ergm.bridge.llr
#'
#' @description `ergm.bridge.dindstart.llk` is a wrapper that uses a
#'   dyad-independent ERGM as a starting point for bridge sampling to
#'   estimate the log-likelihood for a given dyad-dependent model and
#'   parameter configuration.  Note that it only handles binary ERGMs
#'   (`response=NULL`) and with constraints (`constraints=`) that that
#'   do not induce dyadic dependence.
#'
#' @param dind A one-sided formula with the dyad-independent model to use as a
#' starting point. Defaults to the dyad-independent terms found in the formula
#' \code{object} with an overal density term (\code{edges}) added if not
#' redundant.
#' @param coef.dind Parameter configuration for the dyad-independent starting
#' point. Defaults to the MLE of \code{dind}.
#'
#' @return \code{ergm.bridge.dindstart.llk} result list also includes
#'   an `llk` element, with the log-likelihood itself and an
#'   `llk.dind` element, with the log-likelihood of the nearest
#'   dyad-independent model.
#' 
#' @export
ergm.bridge.dindstart.llk<-function(object, response=NULL, constraints=~., coef, obs.constraints=~.-observed, target.stats=NULL, dind=NULL, coef.dind=NULL,  basis=ergm.getnetwork(object), ..., llkonly=TRUE, control=control.ergm.bridge(), verbose=FALSE){
  check.control.class("ergm.bridge", "ergm.bridge.dindstart.llk")
  handle.control.toplevel("ergm.bridge", ...)

  if(verbose) message("Initializing model to obtain the list of dyad-independent terms...")
  ## Here, we need to get the model object to get the list of
  ## dyad-independent terms.
  nw <- basis
  ergm_preprocess_response(nw, response)

  if(is.valued(nw)) stop("Only binary ERGMs are supported at this time.")
  if(!is.null(dind)) stop("Custom dind scaffolding has been disabled. It may be reenabled in the future.")

  m<-ergm_model(object, nw, term.options=control$term.options)
  m.edges <- ergm_model(~edges, nw, term.options = control$term.options)

  if(!is.null(target.stats)){
    if(nparam(m, canonical=TRUE, offset=FALSE)!=length(target.stats)){
      stop("Incorrect length of the target.stats vector: should be ", nparam(m, canonical=TRUE, offset=FALSE), " but is ",length(target.stats),". Note that offset() terms should *not* get target statistics.")
    }
    target.stats <- .align.target.stats.offset(m, target.stats)
    target.stats[is.na(target.stats) & m$etamap$offsetmap] <- summary(m, nw)[is.na(target.stats) & m$etamap$offsetmap]
  }

  q.pos.full <- c(0,cumsum(nparam(m, canonical=FALSE, byterm=TRUE, offset=TRUE)))
  p.pos.full <- c(0,cumsum(nparam(m, canonical=TRUE, byterm=TRUE, offset=FALSE)))
  rng <- function(x, from, to) if(to>=from) x[from:to]
  
  tmp <- .handle.auto.constraints(nw, constraints, obs.constraints, target.stats); nw <- tmp$nw
  if(!is.dyad.independent(ergm_conlist(tmp$conterms,nw,term.options=control$term.options), ergm_conlist(tmp$conterms.obs,nw,term.options=control$term.options))) stop("Bridge sampling with dyad-independent start does not work with dyad-dependent constraints.")

  # If target.stats are given, then we need between passed network and
  # target stats, if any. It also means that the dyad-independent
  # submodel cannot contain any statistics that target.stats does not,
  # so edges are not added on.
  if(!is.null(target.stats) && !is.null(dind)) stop("Non-default dyad-independent model is not supported when target.stats is passed and passed network's statistics do not match it.")

  ## By default, take dyad-independent terms in the formula, fit a
  ## model with these terms and "edges". Terms that are redundant (NA)
  ## get their coefficients zeroed out below.
  ## FIXME: What to do about dyad-independent curved terms?
  offset.dind <- c()
  if(!is.null(target.stats)) ts.dind <- c()
  dindmap <- logical(0)
  if(is.null(dind)){
    terms.full<-list_rhs.formula(object)[!m$term.skipped] # Ensure that terms to be added to the dyad-independent formula are aligned with terms that had actually made it into the model.
    for(i in seq_along(terms.full))
      if(NVL(m$terms[[i]]$dependence, TRUE)){ # Dyad-dependent: drop.
        terms.full[i] <- list(NULL)
        dindmap <- c(dindmap, rep(FALSE, length(m$terms[[i]]$offset)))
      }else{
        dindmap <- c(dindmap, rep(TRUE, length(m$terms[[i]]$offset)))
        if(!is.null(target.stats)) ts.dind <- c(ts.dind, rng(target.stats, p.pos.full[i]+1, p.pos.full[i+1]))
        offset.dind <- c(offset.dind, coef[(q.pos.full[i]+1):q.pos.full[i+1]][m$terms[[i]]$offset]) # Add offset coefficient where applicable.
      }

    terms.full <- c(terms.full, list(as.name("edges")))
    dindmap <- c(dindmap, TRUE)
    if(!is.null(target.stats)) ts.dind <- as.vector(c(ts.dind, edges = network.edgecount(nw)))

    # Copy environment and LHS if present.
    dind <- append_rhs.formula(object[-length(object)], compact(terms.full))

    if(length(object)==3) dind[[2]] <- object[[2]] else dind <- dind[-2]
  }

  message("Fitting the dyad-independent submodel...")
  if(is.null(coef.dind)){
    ergm.dind<-suppressMessages(suppressWarnings(ergm(dind,basis=nw,estimate="MPLE",constraints=constraints,obs.constraints=obs.constraints,eval.loglik=FALSE,control=control.ergm(drop=FALSE, term.options=control$term.options, MPLE.max.dyad.types=control$MPLE.max.dyad.types), offset.coef = offset.dind)))
    etamap.dind <- ergm.dind$etamap
    stats.dind <- ergm.dind$nw.stats

    eta.dind <- ergm.eta(coef(ergm.dind), ergm.dind$etamap)[!ergm.dind$etamap$offsetmap]
    eta.dind <- ifelse(is.na(eta.dind),0,eta.dind)
    llk.dind <- ergm.dind$mple.lik
  }else{
    mple.dind <- suppressMessages(suppressWarnings(ergmMPLE(dind, output="matrix", constraints=constraints,obs.constraints=obs.constraints, control=control.ergm(drop=FALSE, term.options=control$term.options, MPLE.max.dyad.types=control$MPLE.max.dyad.types))))
    etamap.dind <- attr(ergm.dind, "etamap")
    stats.dind <- summary(dind, basis=nw)

    eta.dind <- ergm.eta(coef.dind, etamap.dind)
    lin.pred <- mple.dind$x %*% eta.dind
    llk.dind <- crossprod(lin.pred, mple.dind$response*mple.dind$weights)-sum(log1p(exp(mple.dind$predictor))*mple.dind$weights)
  }

  # If there are target.stats we need to adjust the log-likelihood in
  # case they are different from those to which the dyad-independent
  # submodel was actually fit:
  # l(theta,ts)-l(theta,ns)=sum(theta*(ts-ns)).
  if(!is.null(target.stats)) llk.dind <- llk.dind + c(crossprod(eta.dind, NVL(c(ts.dind), stats.dind[!etamap.dind$offsetmap]) - stats.dind[!etamap.dind$offsetmap]))

  coef.dind <- numeric(length(dindmap))
  coef.dind[dindmap] <- replace(coef(ergm.dind), is.na(coef(ergm.dind)), 0)
  coef.aug <- c(coef, 0)

  form.aug <- append_rhs.formula(object, list(as.name("edges")))

  ## From this point on, target.stats has NAs corresponding to offset
  ## terms.
  if(!is.null(target.stats)) target.stats <- unname(c(target.stats[!m$etamap$offsetmap], ult(ts.dind)))

  if(verbose){
    message("Dyad-independent submodel MLE has likelihood ", format(llk.dind), " at:")
    message_print(coef.dind)
  }

  # NB: Since the LHS network is almost certainly going to be closer
  # to a draw from the full model's MLE than from the submodel's MLE,
  # bridge from the full model to the submodel and subtract below.
  message("Bridging between the dyad-independent submodel and the full model...")
  br <- ergm.bridge.llr(form.aug, constraints = constraints, obs.constraints = obs.constraints,
                        from = coef.aug, to = coef.dind, basis = basis, target.stats = target.stats,
                        control = control, verbose = verbose)
  message("Bridging finished.")
  
  if (llkonly) llk.dind - br$llr
  else c(br, llk.dind = llk.dind, llk = llk.dind - br$llr)
}
