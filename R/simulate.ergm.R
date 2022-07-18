#  File R/simulate.ergm.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
#========================================================================
# This file contains the following 2 functions for simulating ergms
#           <simulate.ergm>
#           <simulate.formula.ergm>
#========================================================================

#' Draw from the distribution of an Exponential Family Random Graph Model
#' 
#' \code{\link[stats]{simulate}} is used to draw from exponential
#' family random network models.  See \code{\link{ergm}} for more
#' information on these models. 
#'
#' 
#'
#' A sample of networks is randomly drawn from the specified model.  The model
#' is specified by the first argument of the function.  If the first argument
#' is a \code{\link{formula}} then this defines the model.  If the first
#' argument is the output of a call to \code{\link{ergm}} then the model used
#' for that call is the one fit -- and unless \code{coef} is specified, the
#' sample is from the MLE of the parameters.  If neither of those are given as
#' the first argument then a Bernoulli network is generated with the
#' probability of ties defined by \code{prob} or \code{coef}.
#' 
#' Note that the first network is sampled after \code{burnin} steps,
#' and any subsequent networks are sampled each \code{interval} steps
#' after the first.
#' 
#' More information can be found by looking at the documentation of
#' \code{\link{ergm}}.
#' 
#' @param object Either a \code{\link{formula}} or an
#' \code{\link{ergm}} object.  The \code{\link{formula}} should be of the form
#' \code{y ~ <model terms>}, where \code{y} is a network object or a matrix
#' that can be coerced to a \code{\link[network]{network}} object.  For the
#' details on the possible \code{<model terms>}, see \code{\link{ergmTerm}}.
#' To create a \code{\link[network]{network}} object in , use the
#' \code{network()} function, then add nodal attributes to it using the
#' \code{\%v\%} operator if necessary.
#' @param nsim Number of networks to be randomly drawn from the given
#' distribution on the set of all networks, returned by the Metropolis-Hastings
#' algorithm.
#' @template seed
#' 
#' @param coef Vector of parameter values for the model from which the
#'   sample is to be drawn.  If \code{object} is of class \code{ergm},
#'   the default value is the vector of estimated coefficients. Can be
#'   set to `NULL` to bypass, but only if `return.args` below is used.
#' 
#' @template response
#' @template reference
#' @template constraints
#'
#' @param observational Inherit observational constraints rather than model
#'   constraints.
#' 
#' @param monitor A one-sided formula specifying one or more terms
#'   whose value is to be monitored. These terms are appended to the
#'   model, along with a coefficient of 0, so their statistics are
#'   returned. An [`ergm_model`] objectcan be passed as well.
#'
#' @template basis
#' 
#' @param statsonly Logical: If TRUE, return only the network statistics, not
#' the network(s) themselves. Deprecated in favor of `output=`.
#' @param esteq Logical: If TRUE, compute the sample estimating equations of an
#' ERGM: if the model is non-curved, all non-offset statistics are returned
#' either way, but if the model is curved, the score estimating function values
#' (3.1) by Hunter and Handcock (2006) are returned instead.
#' 
#' @param output Normally character, one of `"network"` (default),
#'   `"stats"`, `"edgelist"`, or `"ergm_state"`: determines the output
#'   format. Partial matching is performed.
#'
#'   Alternatively, a function with prototype
#'   `function(ergm_state, chain, iter, ...)` that is
#'   called for each returned network, and its return value, rather
#'   than the network itself, is stored. This can be used to, for
#'   example, store the simulated networks to disk without storing
#'   them in memory or compute network statistics not implemented
#'   using the ERGM API, without having to store the networks
#'   themselves.
#'
#' @param simplify Logical: If `TRUE` the output is "simplified":
#'   sampled networks are returned in a single list, statistics from
#'   multiple parallel chains are stacked, etc.. This makes it
#'   consistent with behavior prior to `ergm` 3.10.
#' 
#' @param sequential Logical: If FALSE, each of the \code{nsim} simulated
#' Markov chains begins at the initial network.  If TRUE, the end of one
#' simulation is used as the start of the next.  Irrelevant when \code{nsim=1}.
#'
#' @templateVar mycontrols [control.simulate.ergm()] or [control.simulate.formula()]
#' @template control2
#' @template verbose
#'
#' @param \dots Further arguments passed to or used by methods.
#' 
#' @param return.args Character; if not `NULL`, the `simulate` method
#'   for that particular class will, instead of proceeding for
#'   simulation, instead return its arguments as a list that can be
#'   passed as a second argument to [do.call()] or a lower-level
#'   function such as [ergm_MCMC_sample()]. This can be useful if, for
#'   example, one wants to run several simulations with varying
#'   coefficients and does not want to reinitialize the model and the
#'   proposal every time. Valid inputs at this time are `"formula"`,
#'   "ergm_model", and one of the `"ergm_state"` classes, for the three
#'   respective stopping points.
#'
#' @param do.sim Logical; a deprecated interface superseded by `return.args`,
#'   that saves the inputs to the next level of the function.
#'
#' @return If \code{output=="stats"} an [`mcmc`] object containing the
#'   simulated network statistics. If \code{control$parallel>0}, an
#'   [`mcmc.list`] object. If `simplify=TRUE` (the default), these
#'   would then be "stacked" and converted to a standard [`matrix`]. A
#'   logical vector indicating whether or not the term had come from
#'   the `monitor=` formula is stored in [attr()]-style attribute
#'   `"monitored"`.
#'
#' Otherwise, a representation of the simulated network is returned,
#' in the form specified by `output`. In addition to a network
#' representation or a list thereof, they have the following
#' \code{\link{attr}}-style attributes: \describe{
#'
#' \item{`formula`}{The \code{\link{formula}} used to generate the
#' sample.}
#'
#' \item{`stats`}{An [`mcmc`] or [`mcmc.list`] object as above.}
#'
#' \item{`control`}{Control parameters used to generate the sample.}
#'
#' \item{`constraints`}{Constraints used to generate the sample.}
#'
#' \item{`reference`}{The reference measure for the sample.}
#' 
#' \item{`monitor`}{The monitoring formula.}
#'
#' \item{`response`}{The edge attribute used as a response.}
#'
#' }
#'
#' The following are the permitted network formats: \describe{
#'
#' \item{`"network"`}{If \code{nsim==1}, an object of class
#' \code{network}.  If \code{nsim>1}, it returns an object of class
#' \code{\link{network.list}} (a list of networks) with the
#' above-listed additional attributes.}
#'
#' \item{`"edgelist"`}{An [`edgelist`] representation of the network,
#' or a list thereof, depending on `nsim`.}
#'
#' \item{`"ergm_state"`}{A semi-internal representation of
#' a network consisting of a [`network`] object emptied of edges, with
#' an attached edgelist matrix, or a list thereof, depending on
#' `nsim`.}
#'
#' }
#'
#' If `simplify==FALSE`, the networks are returned as a nested list,
#' with outer list being the parallel chain (including 1 for no
#' parallelism) and inner list being the samples within that chains
#' (including 1, if one network per chain). If `TRUE`, they are
#' concatenated, and if a total of one network had been simulated, the
#' network itself will be returned.
#'
#' @note The actual [`network`] method for [simulate_formula()] is
#'   actually called `.simulate_formula.network()` and is also
#'   exported as an object. This allows it to be overridden by
#'   extension packages, such as `tergm`, but also accessed directly
#'   when needed.
#'
#' @seealso \code{\link{ergm}}, \code{\link[network]{network}},
#'   [ergm_MCMC_sample()] for a demonstration of `return.args=`.
#' @keywords models
#' @examples
#' \dontshow{
#' options(ergm.eval.loglik=FALSE)
#' }
#' #
#' # Let's draw from a Bernoulli model with 16 nodes
#' # and density 0.5 (i.e., coef = c(0,0))
#' #
#' g.sim <- simulate(network(16) ~ edges + mutual, coef=c(0, 0))
#' #
#' # What are the statistics like?
#' #
#' summary(g.sim ~ edges + mutual)
#' #
#' # Now simulate a network with higher mutuality
#' #
#' g.sim <- simulate(network(16) ~ edges + mutual, coef=c(0,2))
#' #
#' # How do the statistics look?
#' #
#' summary(g.sim ~ edges + mutual)
#' #
#' # Let's draw from a Bernoulli model with 16 nodes
#' # and tie probability 0.1
#' #
#' g.use <- network(16,density=0.1,directed=FALSE)
#' #
#' # Starting from this network let's draw 3 realizations
#' # of a edges and 2-star network
#' #
#' g.sim <- simulate(~edges+kstar(2), nsim=3, coef=c(-1.8,0.03),
#'                basis=g.use, control=control.simulate(
#'                  MCMC.burnin=1000,
#'                  MCMC.interval=100))
#' g.sim
#' summary(g.sim)
#' #
#' # attach the Florentine Marriage data
#' #
#' data(florentine)
#' #
#' # fit an edges and 2-star model using the ergm function
#' #
#' gest <- ergm(flomarriage ~ edges + kstar(2))
#' summary(gest)
#' #
#' # Draw from the fitted model (statistics only), and observe the number
#' # of triangles as well.
#' #
#' g.sim <- simulate(gest, nsim=10, 
#'             monitor=~triangles, output="stats",
#'             control=control.simulate.ergm(MCMC.burnin=1000, MCMC.interval=100))
#' g.sim
#'
#' # Custom output: store the edgecount (computed in R), iteration index, and chain index.
#' output.f <- function(x, iter, chain, ...){
#'   list(nedges = network.edgecount(as.network(x)),
#'        chain = chain, iter = iter)
#' }
#' g.sim <- simulate(gest, nsim=3,
#'             output=output.f, simplify=FALSE,
#'             control=control.simulate.ergm(MCMC.burnin=1000, MCMC.interval=100))
#' unclass(g.sim)
#' @name simulate.ergm
#' @importFrom stats simulate
#' @aliases simulate.formula.ergm
#' @export
simulate.formula_lhs_network <- function(object, nsim=1, seed=NULL, ...){
  simulate_formula(object, nsim=nsim, seed=seed, ..., basis=attr(object, ".Basis"))
}

#' @rdname simulate.ergm
#'
#' @export
simulate_formula <- function(object, ..., basis=eval_lhs.formula(object)) {
  UseMethod("simulate_formula", object=basis)  
}

#' @rdname simulate.ergm
#'
#' @rawNamespace S3method(simulate_formula,network,.simulate_formula.network)
#' @aliases simulate_formula.network
#' @method simulate_formula network
#' @export .simulate_formula.network
.simulate_formula.network <- function(object, nsim=1, seed=NULL,
                               coef, response=NULL, reference=~Bernoulli,
                             constraints=~.,
                             observational=FALSE,
                               monitor=NULL,
                               statsonly=FALSE,
                             esteq=FALSE,
                             output=c("network","stats","edgelist","ergm_state"),
                             simplify=TRUE,
                             sequential=TRUE,
                               control=control.simulate.formula(),
                             verbose=FALSE, ..., basis=ergm.getnetwork(object), do.sim=NULL,
                             return.args = NULL){
  if(!missing(do.sim) && !is.null(do.sim)){
    .Deprecate_once(msg=paste0("Use of ",sQuote("do.sim=")," argument has been deprecated. Use ",sQuote("return.args=")," instead."))
    if(!do.sim) return.args <- "ergm_model"
  }
  if(!is.null(return.args) && is(object, return.args))
    return(c(as.list(environment()), list(...)))

  #' @importFrom statnet.common check.control.class
  check.control.class("simulate.formula", myname="ERGM simulate.formula")
  handle.control.toplevel("simulate.formula", ...)

  if(!missing(statsonly)){
    .Deprecate_once(msg=paste0("Use of ",sQuote("statsonly=")," argument has been deprecated. Use ",sQuote("output='stats'")," instead."))
    output <- if(statsonly) "stats" else "network"
  }

  nw <- basis
  nw <- as.network(nw, populate=FALSE)
  ergm_preprocess_response(nw, response)

  mon.m <- if(!is.null(monitor)) as.ergm_model(monitor, nw, term.options=control$term.options)

  # Construct the proposal; this needs to be done here so that the
  # auxiliary requests could be passed to ergm_model().
  if(is(constraints, "ergm_proposal")) proposal <- constraints
  else{
    if(!is.list(constraints)) constraints <- list(constraints)
    constraints <- rep(constraints, length.out=2)
    # Inherit constraints from nw if needed.
    tmp <- .handle.auto.constraints(nw, constraints[[1]], constraints[[2]], NULL)
    nw <- tmp$nw; conterms <- if(observational) tmp$conterms.obs else tmp$conterms

    if (verbose) message("Initializing unconstrained Metropolis-Hastings proposal: ", appendLF=FALSE)
    proposal <- ergm_proposal(conterms, arguments=if(observational) control$obs.MCMC.prop.args else control$MCMC.prop.args,
                              nw=nw, hints=if(observational) control$obs.MCMC.prop else control$MCMC.prop, weights=if(observational) control$obs.MCMC.prop.weights else control$MCMC.prop.weights, class="c",reference=reference, term.options=control$term.options)
    if (verbose) message(sQuote(paste0(proposal$pkgname,":MH_",proposal$name)),".")
  }
  
  # Construct the model.
  if (verbose) message("Initializing model...")
  m <- ergm_model(object, nw, extra.aux=list(proposal=proposal$auxiliaries), term.options=control$term.options)
  proposal$aux.slots <- m$slots.extra.aux$proposal
  if (verbose) message("Model initialized.")

  # Pass the inputs to the simualte method for ergm_model.
    out <- simulate(m, nsim=nsim, seed=seed,
                    coef=coef,
                    constraints=proposal,
                    monitor=mon.m,
                    basis=nw,
                    esteq=esteq,
                    output=output,
                    simplify=simplify,
                    sequential=sequential,
                    control=control,
                    verbose=verbose,
                    return.args=return.args,
                    ...)
    
    if(statsonly || nsim==1) # Then out is either a matrix or a single
                           # network. Return it.
      return(out)
  
    # If we get this far, statsonly==FALSE and nsim > 1, so out is a
    # network.list. Therefore, set the simulation and monitor formulas,
    # which simulate.ergm_model() doesn't know.
    attributes(out) <- c(attributes(out),
                         list(formula=object, monitor=monitor, constraints=constraints, reference=reference))
    out
}

#' @rdname simulate.ergm
#'
#' @export
simulate_formula.ergm_state <- .simulate_formula.network

#' @rdname simulate.ergm
#'
#' @note [simulate.ergm_model()] is a lower-level interface, providing
#'   a [simulate()] method for [`ergm_model`] class. The `basis`
#'   argument is required; `monitor`, if passed, must be an
#'   [`ergm_model`] as well; and `constraints` can be an
#'   [`ergm_proposal`] object instead.
#' @export
simulate.ergm_model <- function(object, nsim=1, seed=NULL,
                                coef, reference=if(is(constraints, "ergm_proposal")) NULL else trim_env(~Bernoulli),
                                constraints=trim_env(~.),
                                observational=FALSE,
                                monitor=NULL,
                                basis=NULL,
                                esteq=FALSE,
                                output=c("network","stats","edgelist","ergm_state"),
                                simplify=TRUE,
                                sequential=TRUE,
                                control=control.simulate.formula(),
                                verbose=FALSE, ..., do.sim=NULL,
                                return.args = NULL){
  if(!missing(do.sim) && !is.null(do.sim)){
    .Deprecate_once(msg=paste0("Use of ",sQuote("do.sim=")," argument has been deprecated. Use ",sQuote("return.args=")," instead."))
    if(!do.sim) return.args <- "ergm_state"
  }
  if(!is.null(return.args) && is(object, return.args))
    return(c(as.list(environment()), list(...)))

  check.control.class(c("simulate.formula", "simulate.ergm_model"), myname="simulate.ergm_model")
  handle.control.toplevel("simulate.formula", ...)
  
  if(!is.null(monitor) && !is(monitor, "ergm_model")) stop("ergm_model method for simulate() requires monitor= argument of class ergm_model or NULL.")
  if(is.null(basis)) stop("ergm_model method for simulate() requires the basis= argument for the initial state of the simulation.")

  if(is.character(output))
    output <- match.arg(output)
  else{
    output.f <- output
    output <- "function"
  }

  # Backwards-compatibility code:
  if("theta0" %in% names(list(...))){
    warning("Passing the parameter vector as theta0= is deprecated. Use coef= instead.")
    coef<-list(...)$theta0
  }
  
  if(!is.null(seed)) {set.seed(as.integer(seed))}
  
  # define nw as either the basis argument or (if NULL) the LHS of the formula
  nw <- basis
  nw0 <- as.network(nw, populate=FALSE)

  m <- c(object, monitor)
  
  # Just in case the user did not give a coef value, set it to zero.
  # (probably we could just return an error in this case!)
  if(missing(coef)) {
    coef <- c(rep(0, nparam(m)))
    warning("No parameter values given, using Bernouli network.")
  }

  if(!is.null(coef)){
    coef <- c(coef, rep(0, nparam(monitor)))
    if(nparam(m)!=length(coef)) stop("coef has ", length(coef) - nparam(monitor), " elements, while the model requires ",nparam(m) - nparam(monitor)," parameters.")
  }

  if(is(constraints, "ergm_proposal")) proposal <- constraints
  else{
    if(is.ergm_state(nw)) warning(sQuote("simulate.ergm_model()"), " has been passed a network in ", sQuote("ergm_state"), " form but not a pre-initialized proposal. Information about missing dyads may be lost.")
    if(!is.list(constraints)) constraints <- list(constraints)
    constraints <- rep(constraints, length.out=2)
    # Inherit constraints from nw if needed.
    tmp <- .handle.auto.constraints(nw0, constraints[[1]], constraints[[2]], NULL)
    nw0 <- tmp$nw; conterms <- if(observational) tmp$conterms.obs else tmp$conterms

    if (verbose) message("Initializing unconstrained Metropolis-Hastings proposal: ", appendLF=FALSE)
    proposal <- ergm_proposal(conterms, arguments=control$MCMC.prop.args,
                              nw=nw0, hints=control$MCMC.prop, weights=control$MCMC.prop.weights, class="c",reference=reference, term.options=control$term.options)
    if (verbose) message(sQuote(paste0(proposal$pkgname,":MH_",proposal$name)),".")
  }

  if(length(proposal$auxiliaries) && !length(m$slots.extra.aux$proposal))
    stop("The proposal appears to be requesting auxiliaries, but the initialized model does not export any proposal auxiliaries.")
  
  # Create vector of current statistics
  curstats<-summary(m, nw, term.options=control$term.options)
  names(curstats) <- param_names(m, canonical=TRUE)

  state <- ergm_state(nw, model=m, proposal=proposal, stats=curstats)

    o <- simulate(state, nsim=nsim, seed=seed,
                  coef,
                  esteq=esteq,
                  output=if(output=="function") output.f else output,
                  simplify=simplify,
                  sequential=sequential,
                  control=control,
                  verbose=verbose,
                  return.args=return.args,
                  ...)

  mon <- rep(c(FALSE,TRUE), c(nparam(m,canonical=!esteq) - NVL3(monitor, nparam(.,canonical=!esteq), 0), NVL3(monitor, nparam(.,canonical=!esteq), 0)))
  if(output=="stats" || is.null(attr(o, "stats"))) attr(o, "monitored") <- mon
  else attr(attr(o, "stats"), "monitored") <- mon
  o
}

#' @describeIn simulate.ergm a low-level function to simulate from an [`ergm_state`] object.
#' @aliases simulate.ergm_state
#' @export
simulate.ergm_state_full <- function(object, nsim=1, seed=NULL,
                                coef,
                                esteq=FALSE,
                                output=c("network","stats","edgelist","ergm_state"),
                                simplify=TRUE,
                                sequential=TRUE,
                                control=control.simulate.formula(),
                                verbose=FALSE, ..., return.args=NULL){
  if(!is.null(return.args)){
    if(is(object, return.args)) return(c(as.list(environment()), list(...)))
    else stop("return.args= is not NULL yet the code has arrived at the actual simulation stage; this likely means an incorrect value had been passed to return.args=")
  }

  if (any(is.nan(coef) | is.na(coef)))
    stop("Illegal value of coef passed to simulate functions")

  if(is.character(output))
    output <- match.arg(output)
  else{
    output.f <- output
    output <- "function"
  }

  state <- object
  m <- as.ergm_model(state)

  # Explain how many iterations and steps will ensue if verbose==TRUE
  if (verbose) message(paste0("Starting MCMC iterations to generate ", nsim,
                              " network", if (nsim > 1) "s"))

  if (nsim > 1) {
    # Only start the cluster if needed. We don't actually need it
    # within this function, but ergm_MCMC_sample() will find it. It
    # should stop automatically when the function exits.
    ergm.getCluster(control, verbose)
  }

  #########################
  ## Main part of function:

  convert_output <-
    switch(output,
           ergm_state = function(states, ...) states,
           network = function(states, ...) lapply(states, as.network),
           edgelist = function(states, ...) lapply(states, as.edgelist),
           "function" = function(states, chain, ...) mapply(output.f, states, chain = chain, iter = seq_along(states), SIMPLIFY = FALSE)
           )

  if(output != "stats") control$MCMC.save_networks <- TRUE
  if(!sequential) control$MCMC.batch <- 1
  else NVL(control$MCMC.batch) <- 0

  if(control$MCMC.batch == 0 || control$MCMC.batch >= nsim){ # also implies sequential==TRUE
    # In this case, we can make one parallelized run of
    # ergm_MCMC_sample.
    control$MCMC.samplesize <- nsim
    z <- ergm_MCMC_sample(state, control, theta=coef, verbose=max(verbose-1,0))
    stats <- z$stats
    if(output != "stats") # then store the returned network:
      nw.list <- mapply(convert_output, z$sampnetworks, chain=seq_along(z$sampnetworks), SIMPLIFY=FALSE)
  }else{
    # Create objects to store output
    if (output!="stats") {
      nw.list <- rep(list(list()),nthreads(control))
    }
    stats <- rep(list(matrix(nrow=0, ncol=nparam(state,canonical=TRUE), 
                             dimnames = list(NULL, param_names(m,canonical=TRUE)))),nthreads(control))
    
    # Call ergm_MCMC_sample once for each batch desired.  This is
    # much slower than when sequential==TRUE and MCMC.batch==0, but
    # here we have a more complicated situation: Either we want a
    # multiple runs (MCMC.batch!=0) or we want to
    # restart each chain at the original network (sequential=FALSE).

    batch_size <- max(control$MCMC.batch, nthreads(control))

    for(i in seq_len(ceiling(nsim/batch_size))){

      control.parallel <- modifyList(control,
                                     list(MCMC.samplesize = batch_size,
                                          MCMC.burnin = if(i==1 || sequential==FALSE) control$MCMC.burnin else control$MCMC.interval))
      z <- ergm_MCMC_sample(state, control.parallel, theta=coef, verbose=max(verbose-1,0))
      
      stats <- mapply(function(s, os) rbind(s, os), stats, z$stats, SIMPLIFY=FALSE)
      
      if(output != "stats"){ # then store the returned network:
        # Concatenate the following...
        newnw <- mapply(convert_output, z$sampnetworks, chain=seq_along(z$sampnetworks), SIMPLIFY=FALSE)
        nw.list <- mapply(c, nw.list, newnw, SIMPLIFY=FALSE)

        if(sequential){ # then update the network state:
          state <- z$networks
        }
      }

      if (verbose) message(sprintf("Finished simulation %d of %d.", i * control$MCMC.batch, nsim))
    }
  }

  if(esteq) stats <- lapply(ergm.estfun, stats, coef, m)

  stats <- as.mcmc.list(lapply(stats, mcmc, start=control$MCMC.burnin+1, thin=control$MCMC.interval))
  if(simplify)
    stats <- as.matrix(stats)[seq_len(nsim),,drop=FALSE]

  if(output=="stats")
    return(stats)
  
  # If we get here, output!="stats".
  if(simplify){
    nw.list <- unlist(nw.list, recursive=FALSE)[seq_len(nsim)] # Concatenate.
  }

  if(length(nw.list)==1&&simplify){
    nw.list <- nw.list[[1]] # Just one network.
  }else{
    attributes(nw.list) <- list(coefficients=coef,
                                control=control,
                                response=names(as.edgelist(
                                  if(is.ergm_state(state)) state
                                  else state[[1]]))[3])
    
    class(nw.list) <- "network.list"
  }
  attr(nw.list, "stats") <- stats

  nw.list
}

#' @rdname simulate.ergm
#'
#' @description The method for [`ergm`] objects inherits the model,
#'   the coefficients, the response attribute, the reference, the
#'   constraints, and most simulation parameters from the model fit,
#'   unless overridden by passing them explicitly. Unless overridden,
#'   the simulation is initialized with a random draw from the fitted
#'   model, saved by [ergm()].
#' 
#' @export
simulate.ergm <- function(object, nsim=1, seed=NULL, 
                          coef=coefficients(object),
                          response=object$network%ergmlhs%"response",
                          reference=object$reference,
                          constraints=list(object$constraints, object$obs.constraints),
                          observational=FALSE,
                          monitor=NULL,
                          basis=object$newnetwork,
                          statsonly=FALSE,
                          esteq=FALSE,
                          output=c("network","stats","edgelist","ergm_state"),
                          simplify=TRUE,
                          sequential=TRUE,
                          control=control.simulate.ergm(),
                          verbose=FALSE, ...) {
  check.control.class(c("simulate.ergm","simulate.formula"), "simulate.ergm")
  handle.control.toplevel("simulate.ergm", ...)

  ### TODO: Figure out when adaptive MCMC controls should be inherited.
  ## control.transfer <-
  ##   c(
  ##     STATIC_MCMC_CONTROLS,
  ##     if(!is.null(control$MCMC.effectiveSize)){
  ##       if(statsonly && sequential && control$MCMC.effectiveSize*1.5<=nsim) ADAPTIVE_MCMC_CONTROLS
  ##       else warning("Adaptive MCMC is only supported when sequential==TRUE, statsonly==TRUE and nsim<=1.5 * target effective size. ", "Adaptive MCMC parameters will be ignored.")
  ##     }
  ##   )

  # If both the passed control and the object's control are NULL (such as if MPLE was estimated), overwrite with simulate.formula()'s defaults.
  formula.control <- control.simulate.formula()
  for(arg in STATIC_MCMC_CONTROLS)
    if(is.null(control[[arg]]))
      control[arg] <- list(NVL(object$control[[arg]], formula.control[[arg]]))

  for(arg in SCALABLE_MCMC_CONTROLS)
    if(is.null(control[[arg]]))
      control[arg] <- list(EVL(object$control[[arg]]*control$MCMC.scale, formula.control[[arg]]))

  control <- set.control.class("control.simulate.formula")
  
  if(!missing(statsonly)){
    .Deprecate_once(msg=paste0("Use of ",sQuote("statsonly=")," argument has been deprecated. Use ",sQuote("output='stats'")," instead."))
    output <- if(statsonly) "stats" else "network"
  }

  simulate(object$formula, nsim=nsim, coef=coef, response=response, reference=reference,
                   esteq=esteq,
                   sequential=sequential, constraints=constraints, observational=observational,
                   monitor=monitor,
                   basis=basis,
           output=output, simplify=simplify,
                   control=control, verbose=verbose, seed=seed, ...)
}
