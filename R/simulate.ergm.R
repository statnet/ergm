#  File R/simulate.ergm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2018 Statnet Commons
#######################################################################
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
#' details on the possible \code{<model terms>}, see \code{\link{ergm-terms}}.
#' To create a \code{\link[network]{network}} object in , use the
#' \code{network()} function, then add nodal attributes to it using the
#' \code{\%v\%} operator if necessary.
#' @param nsim Number of networks to be randomly drawn from the given
#' distribution on the set of all networks, returned by the Metropolis-Hastings
#' algorithm.
#' @template seed
#' @param coef Vector of parameter values for the model from which the sample
#' is to be drawn.  If \code{object} is of class \code{ergm}, the default value
#' is the vector of estimated coefficients.
#' @template response
#' @param reference A one-sided formula specifying the
#' reference measure (\eqn{h(y)}) to be used. (Defaults to \code{~Bernoulli}.)
#' See help for [ERGM reference measures][ergm-references] implemented in
#' the \code{\link[=ergm-package]{ergm}} package.
#' @param constraints A one-sided formula specifying one or more constraints on
#' the support of the distribution of the networks being simulated. See the
#' documentation for a similar argument for \code{\link{ergm}} and see
#' [list of implemented constraints][ergm-constraints] for more information. For
#' \code{simulate.formula}, defaults to no constraints. For
#' \code{simulate.ergm}, defaults to using the same constraints as those with
#' which \code{object} was fitted.
#' 
#' @param monitor A one-sided formula specifying one or more terms
#'   whose value is to be monitored. These terms are appeneded to the
#'   model, along with a coefficient of 0, so their statistics are
#'   returned. An [`ergm_model`] objectcan be passed as well.
#'
#' @param basis An optional \code{\link[network]{network}} object to start the
#' Markov chain.  If omitted, the default is the left-hand-side of the
#' \code{formula}.  If neither a left-hand-side nor a \code{basis} is present,
#' an error results because the characteristics of the network (e.g., size and
#' directedness) must be specified.
#' @param statsonly Logical: If TRUE, return only the network statistics, not
#' the network(s) themselves.
#' @param esteq Logical: If TRUE, compute the sample estimating equations of an
#' ERGM: if the model is non-curved, all non-offset statistics are returned
#' either way, but if the model is curved, the score estimating function values
#' (3.1) by Hunter and Handcock (2006) are returned instead.
#' @param sequential Logical: If FALSE, each of the \code{nsim} simulated
#' Markov chains begins at the initial network.  If TRUE, the end of one
#' simulation is used as the start of the next.  Irrelevant when \code{nsim=1}.
#' @param control A list of control parameters for algorithm tuning.
#' Constructed using \code{\link{control.simulate.ergm}} or
#' \code{\link{control.simulate.formula}}, which have different defaults.
#' @param verbose Logical: If TRUE, extra information is printed as the Markov
#' chain progresses.
#' @param \dots Further arguments passed to or used by methods.
#' @return If \code{statsonly==TRUE} a matrix containing the simulated network
#' statistics. If \code{control$parallel>0}, the statistics from each Markov
#' chain are stacked.
#' 
#' Otherwise, if \code{nsim==1}, an object of class \code{network}.  If
#' \code{nsim>1}, it returns an object of class \code{\link{network.list}}: a
#' list of networks with the following \code{\link{attr}}-style attributes on
#' the list: \item{formula}{The \code{\link{formula}} used to generate the
#' sample.} \item{stats}{The \eqn{\code{nsim}\times p} matrix of network
#' statistics, where \eqn{p} is the number of network statistics specified in
#' the model.} \item{control}{Control parameters used to generate the sample.}
#' \item{constraints}{Constraints used to generate the sample.}
#' \item{reference}{The reference measure for the sample.} \item{monitor}{The
#' monitoring formula.} \item{response}{The edge attribute used as a response.}
#' 
#' If \code{statsonly==FALSE && control$parallel>0} the returned networks are
#' "interleaved", in the sense that for \code{y[i,j]} is the \code{j}th network
#' from MCMC chain \code{i}, the sequence returned if
#' \code{control$parallel==2} is \code{list(y[1,1], y[2,1], y[1,2], y[2,2],
#' y[1,3], y[2,3], ...)}. This is different from the behavior when
#' \code{statsonly==TRUE}. This detail may change in the future.
#' 
#' This object has summary and print methods.
#' @seealso \code{\link{ergm}}, \code{\link[network]{network}}
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
#' # Draw from the fitted model (satatistics only), and observe the number
#' # of triangles as well.
#' #
#' g.sim <- simulate(gest, nsim=10, 
#'             monitor=~triangles, statsonly=TRUE,
#'             control=control.simulate.ergm(MCMC.burnin=1000, MCMC.interval=100))
#' g.sim
#' @name simulate.ergm
#' @importFrom stats simulate
#' @aliases simulate.formula.ergm
#' @S3method simulate formula
#' @export simulate.formula
simulate.formula <- function(object, nsim=1, seed=NULL,
                               coef, response=NULL, reference=~Bernoulli,
                               constraints=~.,
                               monitor=NULL,
                               basis=NULL,
                               statsonly=FALSE,
                               esteq=FALSE,
                               sequential=TRUE,
                               control=control.simulate.formula(),
                             verbose=FALSE, ...) {
  .dep_method("simulate","formula")
  
  #' @importFrom statnet.common check.control.class
  check.control.class("simulate.formula", myname="ERGM simulate.formula")
  control.toplevel(...)
  
  if(!is.null(seed)) {set.seed(as.integer(seed))}
  
  # define nw as either the basis argument or (if NULL) the LHS of the formula
  if (is.null(nw <- basis)) {
    nw <- ergm.getnetwork(object)    
  }
  
  # Do some error-checking on the nw object
  nw <- as.network(nw)
  if(!is.network(nw)){
    stop("A network object on the LHS of the formula or via",
         " the 'basis' argument must be given")
  }
  
  mon.m <- if(!is.null(monitor)) as.ergm_model(monitor, nw, response=response, term.options=control$term.options)

  # Prepare inputs to ergm.getMCMCsample
  m <- c(ergm_model(object, nw, response=response, role="static", term.options=control$term.options), mon.m)
  # Just in case the user did not give a coef value, set it to zero.
  # (probably we could just return an error in this case!)
  if(missing(coef)) {
    coef <- c(rep(0, nparam(m)))
    warning("No parameter values given, using Bernouli network\n\t")
  }

  coef <- c(coef, rep(0, nparam(mon.m)))
  
  if(nparam(m)!=length(coef)) stop("coef has ", length(coef) - nparam(mon.m), " elements, while the model requires ",nparam(m) - nparam(mon.m)," parameters.")

  proposal <- ergm_proposal(constraints,arguments=control$MCMC.prop.args,
                           nw=nw, weights=control$MCMC.prop.weights, class="c",reference=reference,response=response)  

  if (any(is.nan(coef) | is.na(coef)))
    stop("Illegal value of coef passed to simulate.formula")
  
  # Create vector of current statistics
  curstats<-summary(m, nw, response=response, term.options=control$term.options)
  names(curstats) <- param_names(m, canonical=TRUE)

  # prepare control object
  control$MCMC.init.maxedges <- 1+max(control$MCMC.init.maxedges, network.edgecount(nw))
  
  # Explain how many iterations and steps will ensue if verbose==TRUE
  if (verbose) {
    message(paste ("Starting MCMC iterations to generate ", nsim,
                " network", ifelse(nsim>1,"s","")))
  }

  nthreads <- max(
    if(inherits(control$parallel,"cluster")) nrow(summary(control$parallel))
    else control$parallel,
    1)
  
  
  #########################
  ## Main part of function:
  if(sequential && statsonly){
    # In this case, we can make one, parallelized run of
    # ergm.getMCMCsample.
    control$MCMC.samplesize <- nsim
    z <- ergm_MCMC_sample(nw, m, proposal, control, theta=coef, verbose=max(verbose-1,0), response=response)
    # Post-processing: Shift each row by observed statistics.
    out.mat <- sweep(as.matrix(z$stats)[seq_len(nsim),,drop=FALSE], 2, curstats, "+")
  }else{
    # Create objects to store output
    if (!statsonly) { 
      nw.list <- list()
    }
    out.mat <- matrix(nrow=0, ncol=length(curstats), 
                      dimnames = list(NULL, param_names(m,canonical=TRUE))) 
    
    # Call ergm.getMCMCsample once for each network desired.  This is much slower
    # than when sequential==TRUE and statsonly==TRUE, but here we have a 
    # more complicated situation:  Either we want a network for each
    # MCMC iteration (statsonly=FALSE) or we want to restart each chain
    # at the original network (sequential=FALSE).
    if(nthreads>1) curstats <- matrix(curstats, nrow=nthreads, ncol=length(curstats), byrow=TRUE)
    
    # start a cluster if we need to run in parallel more than once
    parallel.toplevel <- NULL     # top level reminder to stop cluster
    clus <- NULL
    if (ceiling(nsim/nthreads) > 1) {
    
      if (inherits(control$parallel,"cluster")) {
        clus <- ergm.getCluster(control, verbose)
      } else if(is.numeric(control$parallel) && control$parallel!=0){
        clus <- ergm.getCluster(control, verbose)
        ergm.cluster.started(FALSE)
        parallel.toplevel <- control$parallel
        control$parallel <- clus
      } else {
        clus <- NULL
        ergm.cluster.started(FALSE)
        if (!is.numeric(control$parallel))
          warning("Unrecognized value passed to parallel control parameter.")
      }
    }
    
    for(i in 1:ceiling(nsim/nthreads)){
      
      control$MCMC.samplesize <- nthreads
      control$MCMC.burnin <- if(i==1 || sequential==FALSE) control$MCMC.burnin else control$MCMC.interval
      z <- ergm_MCMC_sample(nw, m, proposal, control, theta=coef, verbose=max(verbose-1,0), response=response)
      
      out.mat <- rbind(out.mat, curstats + as.matrix(z$stats))
      
      if(!statsonly) # then store the returned network:
        nw.list[[length(nw.list)+1]] <- z$networks[[1]]
      
      if(sequential){ # then update the network state:
        nw <- z$networks
        curstats <- curstats + as.matrix(z$stats)
      }

      if(verbose){message(sprintf("Finished simulation %d of %d.",i, nsim))}
    }
    
    # done with parallel cluster
    if (!is.null(clus)) {
      if (!is.null(parallel.toplevel)) {  # if not NULL, then we started the cluster
        ergm.cluster.started(TRUE)
      }
      ergm.stopCluster(clus)
    }
  }

  out.mat <- out.mat[seq_len(nsim),,drop=FALSE]

  if(esteq) out.mat <- ergm.estfun(out.mat, coef, m)
  
  if (statsonly)
    return(out.mat)
  
  # If we get here, statsonly==FALSE.
  if (nsim==1) {
    return(nw.list[[1]])
  } else {
    nw.list <- nw.list[seq_len(nsim)]
    attributes(nw.list) <- list(formula=object, monitor=monitor,
                                stats=out.mat, coef=coef,
                                control=control,
                                constraints=constraints, reference=reference,
                                monitor=monitor, response=response)

    class(nw.list) <- "network.list"
    return(nw.list)
  }
}

#' @rdname simulate.ergm
#'
#' @description The method for [`ergm`] objects inherits the model,
#'   the coefficients, the response attribute, the reference, the
#'   constraints, and most simulation parameters from the model fit,
#'   unless overridden by passing them explicitly.
#'
#' @note `simulate.ergm()` and `simulate.formula() are currently
#'   exported as functions. This behaviour has been deprecated in
#'   `ergm` 3.9 and will be removed in a future version. Simply use
#'   `simulate()` instead, or [getS3method()] if absolutely necessary.
#'
#' @S3method simulate ergm
#' @export simulate.ergm
simulate.ergm <- function(object, nsim=1, seed=NULL, 
                          coef=object$coef,
                          response=object$response,
                          reference=object$reference,
                          constraints=object$constraints,
                          monitor=NULL,
                          statsonly=FALSE,
                          esteq=FALSE,
                          sequential=TRUE,
                          control=control.simulate.ergm(),
                          verbose=FALSE, ...) {
  .dep_method("simulate","ergm")
  
  check.control.class(c("simulate.ergm","simulate.formula"), "simulate.ergm")
  control.toplevel(...)
  control.transfer <- c("MCMC.burnin", "MCMC.interval", "MCMC.prop.weights", "MCMC.prop.args", "MCMC.packagenames", "MCMC.init.maxedges","parallel","parallel.type","parallel.version.check","term.options")
  for(arg in control.transfer)
    if(is.null(control[[arg]]))
      control[arg] <- list(object$control[[arg]])

  control <- set.control.class("control.simulate.formula")
  
  simulate.formula(object$formula, nsim=nsim, coef=coef, response=response, reference=reference,
                   statsonly=statsonly,
                   esteq=esteq,
                   sequential=sequential, constraints=constraints,
                   monitor=monitor,
                   control=control, verbose=verbose, seed=seed, ...)
}


