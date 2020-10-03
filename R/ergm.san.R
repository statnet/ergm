#  File R/ergm.san.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################

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
#'   [`network`] object.  For the details on the
#'   possible \code{<model terms>}, see \code{\link{ergm-terms}}.  To
#'   create a \code{\link[network]{network}} object in , use the
#'   \code{network()} function, then add nodal attributes to it using
#'   the \code{\%v\%} operator if necessary.
#' @return A network or list of networks that hopefully have network
#'   statistics close to the \code{target.stats} vector.
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
#'   `"edgelist"`, or `"pending_update_network"`: determines the
#'   output format. Partial matching is performed.
#' @param only.last if `TRUE`, only return the last network generated;
#'   otherwise, return a [`network.list`] with `nsim` networks.
#'
#' @param control A list of control parameters for algorithm tuning; see
#' \code{\link{control.san}}.
#' @param verbose Logical or numeric giving the level of verbosity. Higher values produce more verbose output.
#' @param offset.coef A vector of offset coefficients; these must be passed in by the user.  
#' Note that these should be the same set of coefficients one would pass to \code{ergm} via 
#' its \code{offset.coef} argument.  For example, if one has a curved \code{offset(gwesp)} 
#' term in the model, then the \code{decay} parameter should \emph{not} be passed in via 
#' \code{offset.coef}.
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
                        output=c("network","edgelist","pending_update_network"),
                        only.last=TRUE,
                        control=control.san(),
                        verbose=FALSE, 
                        offset.coef=NULL,
                        ...) {
  check.control.class("san", "san")
  control.toplevel(...,myname="san")

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
  # nw is now a network/pending_update_network hybrid class. As long
  # as its edges are only accessed through methods that
  # pending_update_network methods overload, it should be fine.

  proposal<-ergm_proposal(constraints,arguments=control$SAN.prop.args,nw=nw,weights=control$SAN.prop.weights, class="c",reference=reference,response=response)
  model <- ergm_model(formula, nw, response=response, term.options=control$term.options)
    
  ## need to extract post-reviseinit theta offsets here
  model.initial <- ergm_model(formula, nw, response=response, initialfit=TRUE, term.options=control$term.options)  
  
  if(length(offset.coef) != sum(model.initial$etamap$offsettheta)) {
    stop("Length of ", sQuote("offset.coef"), " in SAN is ", length(offset.coef), ", while the number of offset coefficients in the model is ", sum(model.initial$etamap$offsettheta), ".")  
  }
  
  if(any(is.na(offset.coef))) {
    stop("Missing offset coefficients passed to SAN.")
  }
  
  init.coefs <- numeric(nparam(model.initial, canonical=FALSE))
  names(init.coefs) <- param_names(model.initial)
  init.coefs[model.initial$etamap$offsettheta] <- offset.coef
  coefs <- ergm.reviseinit(model, init.coefs)
  offset.coef <- coefs[model$etamap$offsettheta]
    
  san(model, response=response, reference=reference, constraints=proposal, target.stats=target.stats, nsim=nsim, basis=nw, output=output, only.last=only.last, control=control, verbose=verbose, offset.coef=offset.coef, ...)
}

#' @describeIn san A lower-level function that expects a pre-initialized [`ergm_model`].
#' @export
san.ergm_model <- function(object, response=NULL, reference=~Bernoulli, constraints=~., target.stats=NULL,
                           nsim=NULL, basis=NULL,
                           output=c("network","edgelist","pending_update_network"),
                           only.last=TRUE,
                           control=control.san(),
                           verbose=FALSE, 
                           offset.coef=NULL,
                           ...) {
  check.control.class("san", "san")
  control.toplevel(...,myname="san")

  output <- match.arg(output)
  model <- object

  out.list <- list()
  out.mat <- numeric(0)

  if(!is.null(nsim)){
    .Deprecate_once(msg = "nsim= argument for the san() functions has been deprecated. Just use replicate().")
    if(nsim>1 && !is.null(control$seed)) warn("Setting the random seed with nsim>1 will produce a list of identical networks.")
    if(nsim>1){
      return(structure(replicate(nsim,
                                 san(object, response=response, reference=reference, constraints=constraints, target.stats=target.stats,
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
  # nw is now a network/pending_update_network hybrid class. As long
  # as its edges are only accessed through methods that
  # pending_update_network methods overload, it should be fine.

  if(is.null(target.stats)){
    stop("You need to specify target statistic via",
         " the 'target.stats' argument")
  }

  proposal <- if(inherits(constraints, "ergm_proposal")) constraints
              else ergm_proposal(constraints,arguments=control$SAN.prop.args,
                                 nw=nw, weights=control$SAN.prop.weights, class="c",reference=reference,response=response)

  ## we have the post-reviseinit offset thetas
  ## need to remap them to etas using ergm.eta
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
    
  Clist <- ergm.Cprepare(nw, model, response=response)
  
  if (verbose) {
    message(paste("Starting ",control$SAN.maxit," MCMC iteration", ifelse(control$SAN.maxit>1,"s",""),
        " of ", control$SAN.nsteps,
        " steps", ifelse(control$SAN.maxit>1, " each", ""), ".", sep=""))
  }
  maxedges <- max(control$SAN.init.maxedges, Clist$nedges)
  netsumm <- summary(model,nw,response=response)[!offset.indicators]
  target.stats <- vector.namesmatch(target.stats, names(netsumm))
  stats <- netsumm-target.stats
  control$invcov <- diag(1/(nparam(model, canonical=TRUE) - noffset), nparam(model, canonical=TRUE) - noffset)

  nstepss <-
    (if(is.function(control$SAN.nsteps.alloc)) control$SAN.nsteps.alloc(control$SAN.maxit) else control$SAN.nsteps.alloc) %>%
    rep(length.out=control$SAN.maxit) %>%
    (function(x) x/sum(x) * control$SAN.nsteps) %>%
    round()

  z <- NULL
  for(i in 1:control$SAN.maxit){
    if (verbose) {
      message(paste("#", i, " of ", control$SAN.maxit, ": ", sep=""),appendLF=FALSE)
    }
    
    tau <- control$SAN.tau * (if(control$SAN.maxit>1) (1/i-1/control$SAN.maxit)/(1-1/control$SAN.maxit) else 0)
    nsteps <- nstepss[i]
    
    # if we have (essentially) zero temperature, need to zero out the finite offsets for the C code to work properly,
    # even if control$SAN.ignore.finite.offsets is FALSE
    if(abs(tau) < .Machine$double.eps) offsets[is.finite(offsets)] <- 0
    
    z <- ergm_SAN_slave(Clist, proposal, stats, tau, control, verbose,..., prev.run=z, nsteps=nsteps, statindices=statindices, offsetindices=offsetindices, offsets=offsets)

    if(z$status!=0) stop("Error in SAN.")
    
    outnw <- pending_update_network(nw,z,response=response)
    nw <-  outnw
    stats <- z$s[nrow(z$s),]
    # Use *proposal* distribution of statistics for weights.
    invcov <-
      if(control$SAN.invcov.diag) ginv(diag(diag(cov(z$s.prop)), ncol(z$s.prop)), tol=.Machine$double.eps)
      else ginv(cov(z$s.prop), tol=.Machine$double.eps)

    # On resuming, don't accumulate proposal record.
    z$s.prop <- NULL
    # Ensure no statistic has weight 0:
    diag(invcov)[abs(diag(invcov))<.Machine$double.eps] <- min(diag(invcov)[abs(diag(invcov))>=.Machine$double.eps],1)
    invcov <- invcov / sum(diag(invcov)) # Rescale for consistency.
    control$invcov <- invcov
    
    if(verbose){
      message("SAN summary statistics:")
      message_print(target.stats+stats)
      message("Meanstats Goal:")
      message_print(target.stats)
      message("Difference: SAN target.stats - Goal target.stats =")
      message_print(stats)
      message("New statistics scaling =")
      message_print(diag(control$invcov))
      message("Scaled Mahalanobis distance = ", mahalanobis(stats, 0, invcov, inverted=TRUE))
    }
    
    if(!only.last){
      out.list[[i]] <- switch(output,
                              pending_update_network=outnw,
                              network=as.network(outnw),
                              edgelist=as.edgelist(outnw)
                              )
      out.mat <- z$s
      attr(out.mat, "W") <- invcov
    }else{
      if(i<control$SAN.maxit && isTRUE(all.equal(unname(stats), numeric(length(stats))))){
        if(verbose) message("Target statistics matched exactly.")
        break
      }
    }
  }
  if(control$SAN.maxit > 1 && !only.last){
    structure(out.list, formula = formula, networks = out.list, 
              stats = out.mat, class="network.list")
  }else{
    switch(output,
           pending_update_network=outnw,
           network=as.network(outnw),
           edgelist=as.edgelist(outnw)
           )    
  }
}

#' @describeIn ergm-deprecated The developers are not aware of a use case for this function. Please contact them if you would like to prevent its removal.
#' @export
san.ergm <- function(object, formula=object$formula, 
                     constraints=object$constraints, 
                     target.stats=object$target.stats,
                     nsim=NULL, basis=NULL,
                     output=c("network","edgelist","pending_update_network"),
                     only.last=TRUE,
                     control=object$control$SAN.control,
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
              constraints=constraints,
              control=control,
              verbose=verbose, 
              offset.coef=offset.coef, ...)
}

ergm_SAN_slave <- function(Clist,proposal,stats,tau,control,verbose,...,prev.run=NULL, nsteps=NULL, samplesize=NULL, maxedges=NULL, statindices=NULL, offsetindices=NULL, offsets=NULL) {
  if(is.null(prev.run)){ # Start from Clist
    nedges <- c(Clist$nedges,0,0)
    tails <- Clist$tails
    heads <- Clist$heads
    weights <- Clist$weights
  }else{ # Pick up where we left off
    nedges <- prev.run$newnwtails[1]
    tails <- prev.run$newnwtails[2:(nedges+1)]
    heads <- prev.run$newnwheads[2:(nedges+1)]
    weights <- prev.run$newnwweights[2:(nedges+1)]
    nedges <- c(nedges,0,0)
    stats <- prev.run$s[nrow(prev.run$s),]
  }

  if(is.null(nsteps)) nsteps <- control$SAN.nsteps
  if(is.null(samplesize)) samplesize <- control$SAN.samplesize
  if(is.null(maxedges)) maxedges <- control$SAN.init.maxedges

  repeat{
    if(is.null(Clist$weights)){
      z <- .C("SAN_wrapper",
              as.integer(nedges),
              as.integer(tails), as.integer(heads),
              as.integer(Clist$n),
              as.integer(Clist$dir), as.integer(Clist$bipartite),
              as.integer(Clist$nterms),
              as.character(Clist$fnamestring),
              as.character(Clist$snamestring),
              as.character(proposal$name), as.character(proposal$pkgname),
              as.double(c(Clist$inputs,proposal$inputs)), as.double(deInf(tau)),
              s = as.double(c(stats, numeric(length(stats)*(samplesize-1)))), s.prop = double(length(stats)*samplesize),
              as.integer(samplesize),
              as.integer(nsteps), 
              newnwtails = integer(maxedges),
              newnwheads = integer(maxedges),
              as.double(control$invcov),
              as.integer(verbose), as.integer(proposal$arguments$constraints$bd$attribs),
              as.integer(proposal$arguments$constraints$bd$maxout), as.integer(proposal$arguments$constraints$bd$maxin),
              as.integer(proposal$arguments$constraints$bd$minout), as.integer(proposal$arguments$constraints$bd$minin),
              as.integer(proposal$arguments$constraints$bd$condAllDegExact), as.integer(length(proposal$arguments$constraints$bd$attribs)),
              as.integer(maxedges),
              status = integer(1),
              nstats=as.integer(length(statindices)),
              statindices=as.integer(statindices - 1),
              noffsets=as.integer(length(offsetindices)),
              offsetindices=as.integer(offsetindices - 1),
              offsets=as.double(deInf(offsets)),
              PACKAGE="ergm")
      
        # save the results (note that if prev.run is NULL, c(NULL$s,z$s)==z$s.
      z<-list(s=matrix(z$s, ncol=length(statindices), byrow = TRUE), s.prop=matrix(z$s.prop, ncol=length(statindices), byrow = TRUE),
                newnwtails=z$newnwtails, newnwheads=z$newnwheads, status=z$status, maxedges=maxedges)
    }else{
      z <- .C("WtSAN_wrapper",
              as.integer(nedges),
              as.integer(tails), as.integer(heads), as.double(weights),
              as.integer(Clist$n),
              as.integer(Clist$dir), as.integer(Clist$bipartite),
              as.integer(Clist$nterms),
              as.character(Clist$fnamestring),
              as.character(Clist$snamestring),
              as.character(proposal$name), as.character(proposal$pkgname),
              as.double(c(Clist$inputs,proposal$inputs)), as.double(deInf(tau)),
              s = as.double(c(stats, numeric(length(stats)*(samplesize-1)))), s.prop = double(length(stats)*samplesize),
              as.integer(samplesize),
              as.integer(nsteps), 
              newnwtails = integer(maxedges),
              newnwheads = integer(maxedges),
              newnwweights = double(maxedges),
              as.double(control$invcov),
              as.integer(verbose), 
              as.integer(maxedges),
              status = integer(1),
              nstats=as.integer(length(statindices)),
              statindices=as.integer(statindices - 1),
              noffsets=as.integer(length(offsetindices)),
              offsetindices=as.integer(offsetindices - 1),
              offsets=as.double(deInf(offsets)),              
              PACKAGE="ergm")
      # save the results
      z<-list(s=matrix(z$s, ncol=length(statindices), byrow = TRUE), s.prop=matrix(z$s.prop, ncol=length(statindices), byrow = TRUE),
              newnwtails=z$newnwtails, newnwheads=z$newnwheads, newnwweights=z$newnwweights, status=z$status, maxedges=maxedges)
    }

    
    z$s <- rbind(prev.run$s,z$s)
    z$s.prop <- rbind(prev.run$s.prop,z$s.prop)
    
    if(z$status!=1) return(z) # Handle everything except for MCMC_TOO_MANY_EDGES elsewhere.
    
    # The following is only executed (and the loop continued) if too many edges.
    maxedges <- maxedges * 10
    if(!is.null(control$SAN.max.maxedges)){
      if(maxedges == control$SAN.max.maxedges*10) # True iff the previous maxedges exactly equaled control$SAN.max.maxedges and that was too small.
        return(z) # This will kick the too many edges problem upstairs, so to speak.
      maxedges <- min(maxedges, control$SAN.max.maxedges)
    }
  }
}
