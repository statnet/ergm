#  File R/ergm.san.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################

#' Use Simulated Annealing to attempt to match a network to a vector of mean
#' statistics
#' 
#' This function attempts to find a network or networks whose statistics match
#' those passed in via the \code{target.stats} vector.
#' 
#' @details Acceptance probabilities for proposed toggles are computed as 
#'   we now describe.  There are two contributions: one from targeted
#'   statistics and one from offsets.
#'
#' For the targeted statistics, a matrix of weights \code{W} is determined on 
#'   each \code{san} iteration as follows.  On the first iteration, the matrix
#'   \code{W} is the \code{n} by \code{n} identity matrix (\code{n} = number of
#'   target statistics), divided by \code{n}.  On subsequent iterations: if 
#'   \code{control$SAN.invcov.diag = FALSE} (the default), then the matrix 
#'   \code{W} is the inverse of the covariance matrix of the targeted 
#'   statistics, divided by the sum of its (the inverse's) diagonal;
#'   if \code{control$SAN.invcov.diag = TRUE}, then \code{W} is the inverse 
#'   of the diagonal (regarded as a matrix) of the covariance matrix of the 
#'   targeted statistics, divided by the sum of its (the inverse's) diagonal.
#'   In either of these two cases, the covariance matrix is computed based on 
#'   proposals (not acceptances) made on the previous iteration, and the 
#'   normalization for \code{W} is such that \code{sum(diag(W)) = 1}.  The 
#'   component of the acceptance probability coming from the targeted statistics
#'   is then computed for a given \code{W} as \code{exp([y.Wy - x.Wx]/T)} where 
#'   \code{T} is the temperature, \code{y} the column vector of differences 
#'   \code{network statistics - target statistics} computed before the current 
#'   proposal is made, \code{x} the column vector of differences 
#'   \code{network statistics - target statistics} computed assuming the current proposal 
#'   is accepted, and \code{.} the dot product.  If \code{control$SAN.maxit > 1},
#'   then on the \code{i}th iteration, the temperature \code{T} takes the value 
#'   \code{control$SAN.tau * (1/i - 1/control$SAN.maxit)/(1 - 1/control$SAN.maxit)};
#'   if \code{control$SAN.maxit = 1}, then the temperature \code{T} takes the 
#'   value \code{0}.  Thus, \code{T} steps down from \code{control$SAN.tau} to
#'   \code{0} and is always \code{0} on the final iteration.
#'
#' Offsets also contribute to the acceptance probability, as follows.  If 
#'   \code{eta} are the canonical offsets and \code{Delta} the corresponding
#'   change statistics for a given proposal, then the offset contribution to
#'   the acceptance probability is simply \code{exp(eta.Delta)} where 
#'   \code{.} denotes the dot product.  By default, finite offsets are ignored,
#'   but this behavior can be changed by setting
#'   \code{control$SAN.ignore.finite.offsets = FALSE}.
#'
#' The overall acceptance probability is the product of the targeted statistics
#'   contribution and the offset contribution (with the product capped at one).
#' 
#' @param object Either a [`formula`] or an [`ergm`] object. The
#'   [`formula`] should be of the form \code{y ~ <model terms>}, where
#'   \code{y} is a network object or a matrix that can be coerced to a
#'   [`network`] object.  For the details on the possible \code{<model
#'   terms>}, see \code{\link{ergm-terms}}.  To create a
#'   \code{\link[network]{network}} object in , use the
#'   \code{network()} function, then add nodal attributes to it using
#'   the \code{\%v\%} operator if necessary.
#' @return A network or list of networks that hopefully have network
#'   statistics close to the \code{target.stats} vector. Additionally,
#'   [attr()]-style attributes `formula` and `stats` are included.
#' @keywords models
#' @aliases san.default
#' @export
san <- function(object, ...){
 UseMethod("san")
}

#' @noRd
#' @export
san.default <- function(object,...)
{
  stop("Either a ergm object or a formula argument must be given")
}
#' @describeIn san Sufficient statistics are specified by a [`formula`].
#' 
#' @template response
#' @template reference
#' @param formula (By default, the \code{formula} is taken from the \code{ergm}
#' object.  If a different \code{formula} object is wanted, specify it here.
#' @param constraints A one-sided formula specifying one or more constraints on
#' the support of the distribution of the networks being simulated. See the
#' documentation for a similar argument for \code{\link{ergm}} and see
#' [list of implemented constraints][ergm-constraints] for more information. For
#' \code{simulate.formula}, defaults to no constraints. For
#' \code{simulate.ergm}, defaults to using the same constraints as those with
#' which \code{object} was fitted.
#' @param target.stats A vector of the same length as the number of non-offset statistics
#' implied by the formula, which is either \code{object} itself in the case of
#' \code{san.formula} or \code{object$formula} in the case of \code{san.ergm}.
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
  nw <- tmp$nw; constraints <- tmp$constraints

  proposal<-ergm_proposal(constraints,arguments=control$SAN.prop.args,nw=nw, hints=control$SAN.prop, weights=control$SAN.prop.weights, class="c",reference=reference, term.options=control$term.options)
  model <- ergm_model(formula, nw, extra.aux=list(proposal=proposal$auxiliaries), term.options=control$term.options)
  proposal$aux.slots <- model$slots.extra.aux$proposal

  
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

  proposal <- if(inherits(constraints, "ergm_proposal")) constraints
              else{
                # Inherit constraints from nw if needed.
                tmp <- .handle.auto.constraints(nw, constraints, NULL, NULL)
                nw <- tmp$nw; constraints <- tmp$constraints
                ergm_proposal(constraints,arguments=control$SAN.prop.args,
                              nw=nw, hints=control$SAN.prop, weights=control$SAN.prop.weights, class="c",reference=reference, term.options=control$term.options)
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
      if(control$SAN.invcov.diag) ginv(diag(diag(cov(sm.prop)), ncol(sm.prop)), tol=.Machine$double.eps)
      else ginv(cov(sm.prop), tol=.Machine$double.eps)

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

#' @describeIn ergm-deprecated The developers are not aware of a use case for this function. Please contact them if you would like to prevent its removal.
#' @export
san.ergm <- function(object, formula=object$formula, 
                     constraints=object$constraints, 
                     target.stats=object$target.stats,
                     nsim=NULL, basis=NULL,
                     output=c("network","edgelist","ergm_state"),
                     only.last=TRUE,
                     control=object$control$SAN,
                     verbose=FALSE, 
                     offset.coef=NULL,
                     ...) {
  .Deprecate_once('The developers are not aware of a use case for this function. Please contact them if you would like to prevent its removal.')
  output <- match.arg(output)
  san.formula(formula, nsim=nsim, 
              target.stats=target.stats,
              basis=basis,
              reference = object$reference,
              output=output,
              only.last=only.last,
              constraints=constraints,
              control=control,
              verbose=verbose, 
              offset.coef=offset.coef, ...)
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
