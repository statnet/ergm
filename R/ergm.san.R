#  File R/ergm.san.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################

#' Generate networks with a given set of network statistics
#' 
#' This function attempts to find a network or networks whose statistics match
#' those passed in via the `target.stats` vector.
#' 
#' @details The following description is an exegesis of section 4 of Krivitsky
#'   et al. (2022).
#' 
#'   Let \eqn{\mathbf{g}}{g} be a vector of target statistics for the
#'   network we wish to construct. That is, we are given an arbitrary network
#'   \eqn{\mathbf{y}^0 \in \mathcal{Y}}{y0 ∈ Y}, and we seek a network
#'   \eqn{\mathbf{y} \in \mathcal{Y}}{y ∈ Y} such that
#'   \eqn{\mathbf{g}(\mathbf{y}) \approx \mathbf{g}}{g(y) ≈ g} -- ideally equality is achieved,
#'   but in practice we may have to settle for a close approximation. The
#'   variant of simulated annealing is as follows.
#'   
#'   The energy function is defined 
#'   
#'   \deqn{E_W (\mathbf{y}) = (\mathbf{g}(\mathbf{y}) - \mathbf{g})^\mathsf{T} W (\mathbf{g}(\mathbf{y}) - \mathbf{g}),}{E_W (y) = (g(y) - g)^T W (g(y) - g),}
#'   
#'   with \eqn{W} a symmetric positive (barring multicollinearity in statistics)
#'   definite matrix of weights. This function achieves 0 only if the target is
#'   reached. A good choice of this matrix yields a more efficient search.
#'
#'   A standard simulated annealing loop is used, as described below, with some
#'   modifications. In particular, we allow the user to specify a vector of
#'   offsets \eqn{\eta}{η} to bias the annealing, with \eqn{\eta_k = 0}{η_k = 0} 
#'   denoting no offset. Offsets can be used with SAN to forbid certain
#'   statistics from ever increasing or decreasing. As with [ergm()], offset
#'   terms are specified using the [offset()] decorator and their coefficients
#'   specified with the `offset.coef` argument. By default, finite offsets are
#'   ignored by, but this can be overridden by setting the [control.san()]
#'   argument `SAN.ignore.finite.offsets = FALSE`.
#'   
#'   The number of simulated annealing runs is specified by the `SAN.maxit`
#'   control parameter and the initial value of the temperature \eqn{T} is set
#'   to `SAN.tau`. The value of \eqn{T} decreases linearly until \eqn{T = 0} 
#'   at the last run, which implies that all proposals that increase 
#'   \eqn{E_W (\mathbf{y})}{E_W(y)} are rejected. The weight matrix \eqn{W} 
#'   is initially set to \eqn{I_p / p}, where \eqn{I_p} is the identity matrix
#'   of an appropriate dimension. For weight \eqn{W} and temperature \eqn{T},
#'   the simulated annealing iteration proceeds as follows:
#'   
#'   1. Test if \eqn{E_W(\mathbf{y}) = 0}{E_W(y) = 0}. If so, then exit.
#'   2. Generate a perturbed network \eqn{\mathbf{y^*}}{y*} from a proposal that
#'      respects the model constraints. (This is typically the same proposal as
#'      that used for MCMC.)
#'   3. Store the quantity 
#'      \eqn{\mathbf{g}(\mathbf{y^*}) - \mathbf{g}(\mathbf{y})}{g(y*) - g(y)} 
#'      for later use.
#'   4. Calculate acceptance probability
#' 
#'      \deqn{\alpha = \exp[ - (E_W (\mathbf{y^*}) - E_W (\mathbf{y})) / T + \eta^\mathsf{T} (\mathbf{g}(\mathbf{y^*}) - \mathbf{g}(\mathbf{y}))]}{α = exp( - E_W(y*) - E_W(y) / T + η' (g(y*) - g(y)) ).}
#'    
#'      (If \eqn{|\eta_k| = \infty}{|η_k| = ∞} and \eqn{g_k (\mathbf{y^*}) - g_k (\mathbf{y}) = 0}{g_k (y) - g_k (y) = 0}, their product is defined to be 0.)
#'   5. Replace \eqn{\mathbf{y}}{y} with \eqn{\mathbf{y^*}}{y} with probability
#'      \eqn{\min(1, \alpha)}{min(1, α)}.
#'
#'
#'   After the specified number of iterations, \eqn{T} is updated as described
#'   above, and \eqn{W} is recalculated by first computing a matrix \eqn{S}, the
#'   sample covariance matrix of the proposed differences stored in Step 3
#'   (i.e., whether or not they were rejected), then
#'   \eqn{W = S^+ / tr(S^+)}{W = S+ / tr(S+)}, where \eqn{S^+}{S+} is the
#'   Moore–Penrose pseudoinverse of \eqn{S} and \eqn{tr(S^+)}{tr(S+)} is the
#'   trace of \eqn{S^+}{S+}. The differences in Step 3 closely reflect the
#'   relative variances and correlations among the network statistics.
#' 
#'   In Step 2, the many options for MCMC proposals can provide for effective
#'   means of speeding the SAN algorithm's search for a viable network.
#' 
#' @param object Either a [`formula`] or some other supported
#'   representation of an ERGM, such as an [`ergm_model`] object.
#'   [`formula`] should be of the form \code{y ~ <model terms>}, where
#'   \code{y} is a network object or a matrix that can be coerced to a
#'   [`network`] object.  For the details on the possible \code{<model
#'   terms>}, see \code{\link{ergmTerm}}.  To create a
#'   \code{\link[network]{network}} object in , use the
#'   \code{network()} function, then add nodal attributes to it using
#'   the \code{\%v\%} operator if necessary.
#' 
#' @return A network or list of networks that hopefully have network
#'   statistics close to the \code{target.stats} vector. No guarantees
#'   are provided about their probability distribution. Additionally,
#'   [attr()]-style attributes `formula` and `stats` are included.
#' 
#' @references Krivitsky, P. N., Hunter, D. R., Morris, M., & Klumb, C. (2022).
#'   ergm 4: Computational Improvements. arXiv preprint arXiv:2203.08198.
#' 
#' 
#' @keywords models
#' @aliases san.default
#' @export
san <- function(object, ...){
 UseMethod("san")
}

#' @describeIn san Sufficient statistics are specified by a [`formula`].
#' 
#' @template response
#' @template reference
#' @param formula (By default, the \code{formula} is taken from the \code{ergm}
#' object.  If a different \code{formula} object is wanted, specify it here.
#' @template constraints
#' @param target.stats A vector of the same length as the number of non-offset statistics
#' implied by the formula.
#' @param nsim Number of networks to generate. Deprecated: just use [replicate()].
#' @param basis If not NULL, a \code{network} object used to start the Markov
#' chain.  If NULL, this is taken to be the network named in the formula.
#'
#' @param output Character, one of `"network"` (default),
#'   `"edgelist"`, or `"ergm_state"`: determines the
#'   output format. Partial matching is performed.
#' @param only.last if `TRUE`, only return the last network generated;
#'   otherwise, return a [`network.list`] with `nsim` networks.
#'
#' @templateVar mycontrol control.san
#' @template control
#' @template verbose
#'
#' @param offset.coef A vector of offset coefficients; these must be passed in by the user.  
#' Note that these should be the same set of coefficients one would pass to \code{ergm} via 
#' its \code{offset.coef} argument.
#' @param \dots Further arguments passed to other functions.
#' @examples
#' \donttest{
#' # initialize x to a random undirected network with 50 nodes and a density of 0.1
#' x <- network(50, density = 0.05, directed = FALSE)
#'  
#' # try to find a network on 50 nodes with 300 edges, 150 triangles,
#' # and 1250 4-cycles, starting from the network x
#' y <- san(x ~ edges + triangles + cycle(4), target.stats = c(300, 150, 1250))
#' 
#' # check results
#' summary(y ~ edges + triangles + cycle(4))
#' 
#' # initialize x to a random directed network with 50 nodes
#' x <- network(50)
#' 
#' # add vertex attributes
#' x %v% 'give' <- runif(50, 0, 1)
#' x %v% 'take' <- runif(50, 0, 1)
#' 
#' # try to find a set of 100 directed edges making the outward sum of
#' # 'give' and the inward sum of 'take' both equal to 62.5, so in
#' # edges (i,j) the node i tends to have above average 'give' and j
#' # tends to have above average 'take'
#' y <- san(x ~ edges + nodeocov('give') + nodeicov('take'), target.stats = c(100, 62.5, 62.5))
#' 
#' # check results
#' summary(y ~ edges + nodeocov('give') + nodeicov('take'))
#' 
#' 
#' # initialize x to a random undirected network with 50 nodes
#' x <- network(50, directed = FALSE)
#' 
#' # add a vertex attribute
#' x %v% 'popularity' <- runif(50, 0, 1)
#' 
#' # try to find a set of 100 edges making the total sum of
#' # popularity(i) and popularity(j) over all edges (i,j) equal to
#' # 125, so nodes with higher popularity are more likely to be
#' # connected to other nodes
#' y <- san(x ~ edges + nodecov('popularity'), target.stats = c(100, 125))
#'  
#' # check results
#' summary(y ~ edges + nodecov('popularity'))
#' 
#' # creates a network with denser "core" spreading out to sparser
#' # "periphery"
#' plot(y)
#' }
#' @export
san.formula <- function(object, response=NULL, reference=~Bernoulli, constraints=~., target.stats=NULL,
                        nsim=NULL, basis=NULL,
                        output=c("network","edgelist","ergm_state"),
                        only.last=TRUE,
                        control=control.san(),
                        verbose=FALSE, 
                        offset.coef=NULL,
                        ...) {
  check.control.class("san", "san")
  handle.control.toplevel("san", ...)

  output <- match.arg(output)

  formula <- object

  if(!is.null(basis)) {
    nw <- basis
  } else {
    nw <- ergm.getnetwork(formula)
  }
  if(inherits(nw,"network.list")){
    nw <- nw$networks[[1]]
  }
  if(is.null(target.stats)){
    stop("You need to specify target statistic via",
         " the 'target.stats' argument")
  }

  nw <- as.network(ensure_network(nw), populate=FALSE)
  # nw is now a network/ergm_state hybrid class. As long
  # as its edges are only accessed through methods that
  # ergm_state methods overload, it should be fine.
  if(!is(nw,"ergm_state")) ergm_preprocess_response(nw, response)

  # Inherit constraints from nw if needed.
  tmp <- .handle.auto.constraints(nw, constraints, NULL, NULL)
  nw <- tmp$nw; conterms <- tmp$conterms

  if (verbose) message("Initializing unconstrained Metropolis-Hastings proposal: ", appendLF=FALSE)
  proposal<-ergm_proposal(conterms,arguments=control$SAN.prop.args,nw=nw, hints=control$SAN.prop, weights=control$SAN.prop.weights, class="c",reference=reference, term.options=control$term.options)
  if (verbose) message(sQuote(paste0(proposal$pkgname,":MH_",proposal$name)),".")
  if (verbose) message("Initializing model...")
  model <- ergm_model(formula, nw, extra.aux=list(proposal=proposal$auxiliaries), term.options=control$term.options)
  proposal$aux.slots <- model$slots.extra.aux$proposal
  if (verbose) message("Model initialized.")

  
  if(length(offset.coef) != sum(model$etamap$offsettheta)) {
    stop("Length of ", sQuote("offset.coef"), " in SAN is ", length(offset.coef), ", while the number of offset coefficients in the model is ", sum(model$etamap$offsettheta), ".")  
  }
  
  if(any(is.na(offset.coef))) {
    stop("Missing offset coefficients passed to SAN.")
  }
  
  san(model, reference=reference, constraints=proposal, target.stats=target.stats, nsim=nsim, basis=nw, output=output, only.last=only.last, control=control, verbose=verbose, offset.coef=offset.coef, ...)
}

#' @describeIn san A lower-level function that expects a pre-initialized [`ergm_model`].
#' @export
san.ergm_model <- function(object, reference=~Bernoulli, constraints=~., target.stats=NULL,
                           nsim=NULL, basis=NULL,
                           output=c("network","edgelist","ergm_state"),
                           only.last=TRUE,
                           control=control.san(),
                           verbose=FALSE, 
                           offset.coef=NULL,
                           ...) {
  check.control.class("san", "san")
  handle.control.toplevel("san", ...)

  output <- match.arg(output)
  model <- object

  out.list <- list()
  out.mat <- numeric(0)

  if(!is.null(nsim)){
    .Deprecate_once(msg = "nsim= argument for the san() functions has been deprecated. Just use replicate().")
    if(nsim>1 && !is.null(control$seed)) warn("Setting the random seed with nsim>1 will produce a list of identical networks.")
    if(nsim>1){
      return(structure(replicate(nsim,
                                 san(object, reference=reference, constraints=constraints, target.stats=target.stats,
                                     basis=basis,
                                     output=output,
                                     only.last=only.last,
                                     control=control,
                                     verbose=verbose,
                                     offset.coef=offset.coef,
                                     ...)),
                       class="network.list"))
    }
  }

  
  if(!is.null(control$seed)) set.seed(as.integer(control$seed))
  nw <- basis
  nw <- as.network(ensure_network(nw), populate=FALSE)
  # nw is now a network/ergm_state hybrid class. As long
  # as its edges are only accessed through methods that
  # ergm_state methods overload, it should be fine.

  if(is.null(target.stats)){
    stop("You need to specify target statistic via",
         " the 'target.stats' argument")
  }

  if(inherits(constraints, "ergm_proposal")) proposal <- constraints
  else{
    # Inherit constraints from nw if needed.
    tmp <- .handle.auto.constraints(nw, constraints, NULL, NULL)
    nw <- tmp$nw; conterms <- tmp$conterms
    if (verbose) message("Initializing unconstrained Metropolis-Hastings proposal: ", appendLF=FALSE)
    proposal <- ergm_proposal(conterms,arguments=control$SAN.prop.args,
                              nw=nw, hints=control$SAN.prop, weights=control$SAN.prop.weights, class="c",reference=reference, term.options=control$term.options)
    if (verbose) message(sQuote(paste0(proposal$pkgname,":MH_",proposal$name)),".")
  }

  if(length(proposal$auxiliaries) && !length(model$slots.extra.aux$proposal))
    stop("The proposal appears to be requesting auxiliaries, but the initialized model does not export any proposal auxiliaries.")

  ## need to remap thetas to etas using ergm.eta
  ## then just keep the offsets because we don't
  ## care about the non-offset coefs
  coefs <- numeric(nparam(model, canonical=FALSE))
  coefs[model$etamap$offsettheta] <- offset.coef
  etas <- ergm.eta(coefs, model$etamap)
  
  
  offset.indicators <- model$etamap$offsetmap
  
  offsetindices <- which(offset.indicators)
  statindices <- which(!offset.indicators)
  offsets <- etas[offset.indicators]
  if(control$SAN.ignore.finite.offsets) offsets[is.finite(offsets)] <- 0
  
  noffset <- sum(offset.indicators)
  
  if (verbose) {
    message(paste("Starting ",control$SAN.maxit," SAN iteration", ifelse(control$SAN.maxit>1,"s",""),
        " of ", control$SAN.nsteps,
        " steps", ifelse(control$SAN.maxit>1, " each", ""), ".", sep=""))
  }
  netsumm<-summary(model,nw)[!offset.indicators]
  target.stats <- vector.namesmatch(target.stats, names(netsumm))
  stats <- netsumm-target.stats
  invcov.dim <- nparam(model, canonical=TRUE) - noffset
  NVL(control$SAN.invcov) <- diag(1/invcov.dim, invcov.dim)
  if(!is.SPD(control$SAN.invcov) || nrow(control$SAN.invcov) != invcov.dim)
    stop(sQuote("control$SAN.invcov"), " parameter should be a square matrix of dimension equal to the number of non-offset statistics")

  nstepss <-
    (if(is.function(control$SAN.nsteps.alloc)) control$SAN.nsteps.alloc(control$SAN.maxit) else control$SAN.nsteps.alloc) %>%
    rep(length.out=control$SAN.maxit) %>%
    (function(x) x/sum(x) * control$SAN.nsteps) %>%
    round()

  state <- ergm_state(nw, model=model, proposal=proposal, stats=stats)
  sm <- NULL
  for(i in 1:control$SAN.maxit){
    if (verbose) {
      message(paste("#", i, " of ", control$SAN.maxit, ": ", sep=""),appendLF=FALSE)
    }
    
    tau <- control$SAN.tau * (if(control$SAN.maxit>1) (1/i-1/control$SAN.maxit)/(1-1/control$SAN.maxit) else 0)
    nsteps <- nstepss[i]
    
    # if we have (essentially) zero temperature, need to zero out the finite offsets for the C code to work properly,
    # even if control$SAN.ignore.finite.offsets is FALSE
    if(abs(tau) < .Machine$double.eps) offsets[is.finite(offsets)] <- 0
    
    z <- ergm_SAN_slave(state, tau, control, verbose,..., nsteps=nsteps, statindices=statindices, offsetindices=offsetindices, offsets=offsets)
    state <- z$state
    sm <- rbind(sm, z$s)
    sm.prop <- z$s.prop

    if(z$status!=0) stop("Error in SAN.")
    
    stats <- sm[nrow(sm),]
    # Use *proposal* distribution of statistics for weights.
    invcov <-
      if(control$SAN.invcov.diag) sginv(diag(diag(cov(sm.prop)), ncol(sm.prop)), tol=.Machine$double.eps)
      else sginv(cov(sm.prop), tol=.Machine$double.eps)

    # Ensure no statistic has weight 0:
    diag(invcov)[abs(diag(invcov))<.Machine$double.eps] <- min(diag(invcov)[abs(diag(invcov))>=.Machine$double.eps],1)
    invcov <- invcov / sum(diag(invcov)) # Rescale for consistency.
    control$SAN.invcov <- invcov
    
    if(verbose){
      message("SAN summary statistics:")
      message_print(target.stats+stats)
      message("Meanstats Goal:")
      message_print(target.stats)
      message("Difference: SAN target.stats - Goal target.stats =")
      message_print(stats)
      message("New statistics scaling =")
      message_print(diag(control$SAN.invcov))
      message("Scaled Mahalanobis distance = ", mahalanobis(stats, 0, invcov, inverted=TRUE))
    }
    
    out.mat <- z$s
    attr(out.mat, "W") <- invcov
    if(!only.last){
      out.list[[i]] <- switch(output,
                              ergm_state=state,
                              network=as.network(state),
                              edgelist=as.edgelist(state)
                              )
    }else{
      if(i<control$SAN.maxit && isTRUE(all.equal(unname(stats), numeric(length(stats))))){
        if(verbose) message("Target statistics matched exactly.")
        break
      }
    }
  }
  if(control$SAN.maxit > 1 && !only.last){
    structure(out.list, formula = formula,
              stats = out.mat, class="network.list")
  }else{
    structure(
      switch(output,
             ergm_state=state,
             network=as.network(state),
             edgelist=as.edgelist(state)
             ),
      stats = out.mat
    )
  }
}

#' Internal Function to Perform Simulated Annealing
#'
#' This is an internal function, not normally called directly by the
#' user. The \code{ergm_SAN_slave} function samples networks and
#' network statistics using a simulated annealing (SAN) algorithm via
#' \code{SAN_wrapper}.
#' 
#' @param state an [`ergm_state`] representing the sampler state, containing information about the network, the model, the proposal, and current statistics.
#'
#' @templateVar mycontrol control.san
#' @param tau a scalar; temperature to use; higher temperature means more proposals that "worsen" the statistics are accepted.
#' @param nsteps an integer; number of SAN proposals.
#' @param samplesize an integer; number of network statistics to return.
#' @param statindices,offsetindices,offsets specification for offset handling; see [san.formula()] implementation.
#' @template control
#' @template verbose
#' @param ... additional arguments, currently unused.
#' @keywords internal
#' @export
ergm_SAN_slave <- function(state, tau,control,verbose, ..., nsteps=NULL, samplesize=NULL, statindices=NULL, offsetindices=NULL, offsets=NULL){
  on.exit(ergm_Cstate_clear())

  state$proposal$flags$SAN <- TRUE
  if(is.null(nsteps)) nsteps <- control$SAN.nsteps
  if(is.null(samplesize)) samplesize <- control$SAN.samplesize

  z <-
    if(!is.valued(state)){
      .Call("SAN_wrapper",
            state,
            # SAN settings
            as.double(deInf(tau)),
            as.integer(samplesize),
            as.integer(nsteps),
            as.double(control$SAN.invcov),
            statindices=as.integer(statindices - 1),
            offsetindices=as.integer(offsetindices - 1),
            offsets=as.double(deInf(offsets)),
            as.integer(verbose),
            PACKAGE="ergm")
    }else{
      .Call("WtSAN_wrapper",
            state,
            # SAN settings
            as.double(deInf(tau)),
            as.integer(samplesize),
            as.integer(nsteps),
            as.double(control$SAN.invcov),
            statindices=as.integer(statindices - 1),
            offsetindices=as.integer(offsetindices - 1),
            offsets=as.double(deInf(offsets)),
            as.integer(verbose),
            PACKAGE="ergm")
    }
  # save the results
  z$s <- matrix(z$s, ncol=length(statindices), byrow = TRUE)
  z$s.prop <- matrix(z$s.prop, ncol=length(statindices), byrow = TRUE)
  colnames(z$s) <- colnames(z$s.prop) <- param_names(state, canonical=TRUE)[statindices]
  z$state <- update(z$state)

  z
}



# Acceptance probabilities for proposed toggles are computed as
#   we now describe.  There are two contributions: one from targeted
#   statistics and one from offsets.
# 
# For the targeted statistics, a matrix of weights \code{W} is determined on
#   each \code{san} iteration as follows.  On the first iteration, the matrix
#   \code{W} is the \code{n} by \code{n} identity matrix (\code{n} = number of
#   target statistics), divided by \code{n}.  On subsequent iterations: if
#   \code{control$SAN.invcov.diag = FALSE} (the default), then the matrix
#   \code{W} is the inverse of the covariance matrix of the targeted
#   statistics, divided by the sum of its (the inverse's) diagonal;
#   if \code{control$SAN.invcov.diag = TRUE}, then \code{W} is the inverse
#   of the diagonal (regarded as a matrix) of the covariance matrix of the
#   targeted statistics, divided by the sum of its (the inverse's) diagonal.
#   In either of these two cases, the covariance matrix is computed based on
#   proposals (not acceptances) made on the previous iteration, and the
#   normalization for \code{W} is such that \code{sum(diag(W)) = 1}.  The
#   component of the acceptance probability coming from the targeted statistics
#   is then computed for a given \code{W} as \code{exp([y.Wy - x.Wx]/T)} where
#   \code{T} is the temperature, \code{y} the column vector of differences
#   \code{network statistics - target statistics} computed before the current
#   proposal is made, \code{x} the column vector of differences
#   \code{network statistics - target statistics} computed assuming the current proposal
#   is accepted, and \code{.} the dot product.  If \code{control$SAN.maxit > 1},
#   then on the \code{i}th iteration, the temperature \code{T} takes the value
#   \code{control$SAN.tau * (1/i - 1/control$SAN.maxit)/(1 - 1/control$SAN.maxit)};
#   if \code{control$SAN.maxit = 1}, then the temperature \code{T} takes the
#   value \code{0}.  Thus, \code{T} steps down from \code{control$SAN.tau} to
#   \code{0} and is always \code{0} on the final iteration.
# 
# Offsets also contribute to the acceptance probability, as follows.  If
#   \code{eta} are the canonical offsets and \code{Delta} the corresponding
#   change statistics for a given proposal, then the offset contribution to
#   the acceptance probability is simply \code{exp(eta.Delta)} where
#   \code{.} denotes the dot product.  By default, finite offsets are ignored,
#   but this behavior can be changed by setting
#   \code{control$SAN.ignore.finite.offsets = FALSE}.
# 
# The overall acceptance probability is the product of the targeted statistics
#   contribution and the offset contribution (with the product capped at one).
