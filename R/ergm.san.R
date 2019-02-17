#  File R/ergm.san.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################

#' Use Simulated Annealing to attempt to match a network to a vector of mean
#' statistics
#' 
#' This function attempts to find a network or networks whose statistics match
#' those passed in via the \code{target.stats} vector.
#' 
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
#' @param response Name of the edge attribute whose value
#' is to be modeled. Defaults to \code{NULL} for simple presence or absence.
#' @param reference One-sided formula whose RHS gives the
#' reference measure to be used. (Defaults to \code{~Bernoulli}.)
#' @param formula (By default, the \code{formula} is taken from the \code{ergm}
#' object.  If a different \code{formula} object is wanted, specify it here.
#' @param constraints A one-sided formula specifying one or more constraints on
#' the support of the distribution of the networks being simulated. See the
#' documentation for a similar argument for \code{\link{ergm}} and see
#' [list of implemented constraints][ergm-constraints] for more information. For
#' \code{simulate.formula}, defaults to no constraints. For
#' \code{simulate.ergm}, defaults to using the same constraints as those with
#' which \code{object} was fitted.
#' @param target.stats A vector of the same length as the number of terms
#' implied by the formula, which is either \code{object} itself in the case of
#' \code{san.formula} or \code{object$formula} in the case of \code{san.ergm}.
#' @param nsim Number of desired networks.
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
#' @param \dots Further arguments passed to other functions.
#' @examples
#' # initialize x to a random undirected network with 100 nodes and a density of 0.1
#' x <- network(100, density = 0.1, directed = F)
#'  
#' # try to find a network on 100 nodes with 600 edges, 300 triangles, and 2500 4-cycles, starting from the network x
#' y <- san(x ~ edges + triangles + cycle(4), target.stats = c(600, 300, 2500))
#' 
#' summary(y ~ edges + triangles + cycle(4))
#'    edges triangle   cycle4 
#'      596      301     2502 
#' # close, but we didn't quite get there after one iteration
#' 
#' # try chaining 5 iterations together, returning 5 networks as a network.list
#' z <- san(x ~ edges + triangles + cycle(4), target.stats = c(600, 300, 2500), nsim = 5, only.last=F)
#' 
#' summary(z ~ edges + triangles + cycle(4))
#'      edges triangle cycle4
#' [1,]   597      301   2501
#' [2,]   597      299   2500
#' [3,]   599      299   2498
#' [4,]   600      300   2500
#' [5,]   600      300   2500
#' # got there after 4 iterations
#' 
#' 
#' # initialize x to a random directed network with 100 nodes
#' x <- network(100)
#' 
#' # add vertex attributes
#' x %v% 'give' <- runif(100, 0, 1)
#' x %v% 'take' <- runif(100, 0, 1)
#' 
#' # try to find a set of 200 directed edges making the outward sum of 'give' and the inward sum of 'take' both equal to 125,
#' # so in edges (i,j) the node i tends to have above average 'give' and j tends to have above average 'take'
#' y <- san(x ~ edges + nodeocov('give') + nodeicov('take'), target.stats = c(200, 125, 125))
#' 
#' summary(y ~ edges + nodeocov('give') + nodeicov('take'))
#'         edges nodeocov.give nodeicov.take 
#'       200.000       124.183       124.529 
#' # pretty close; we can't expect perfection on this one
#' 
#' 
#' # initialize x to a random undirected network with 100 nodes
#' x <- network(100, directed = F)
#' 
#' # add a vertex attribute
#' x %v% 'popularity' <- runif(100, 0, 1)
#' 
#' # try to find a set of 200 edges making the total sum of popularity(i) and popularity(j) over all edges (i,j) equal to 250,
#' # so nodes with higher popularity are more likely to be connected to other nodes
#' y <- san(x ~ edges + nodecov('popularity'), target.stats = c(200, 250))
#'  
#' summary(y ~ edges + nodecov('popularity'))
#'              edges nodecov.popularity 
#'           200.0000           249.4979 
#' 
#' # creates a network with denser "core" spreading out to sparser "periphery"
#' plot(y)
#' @export
san.formula <- function(object, response=NULL, reference=~Bernoulli, constraints=~., target.stats=NULL,
                        nsim=1, basis=NULL,
                        output=c("network","edgelist","pending_update_network"),
                        only.last=TRUE,
                        control=control.san(),
                        verbose=FALSE, ...) {
  check.control.class("san", "san")
  control.toplevel(...,myname="san")

  output <- match.arg(output)

  out.list <- list()
  out.mat <- numeric(0)
  formula <- object

  if(!is.null(control$seed)) set.seed(as.integer(control$seed))
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

  model <- ergm_model(formula, nw, response=response, term.options=control$term.options)
  Clist <- ergm.Cprepare(nw, model, response=response)
  
  proposal<-ergm_proposal(constraints,arguments=control$SAN.prop.args,nw=nw,weights=control$SAN.prop.weights, class="c",reference=reference,response=response)
# if(is.null(control$coef)) {
#   warning("No parameter values given, using the MPLE for the passed network.\n\t")
# }
# control$coef <- c(control$coef[1],rep(0,Clist$nstats-1))
    
  maxedges <- max(control$SAN.init.maxedges, Clist$nedges)
  netsumm<-summary(model,nw,response=response)
  target.stats <- vector.namesmatch(target.stats, names(netsumm))
  stats <- netsumm-target.stats
  control$invcov <- diag(1/nparam(model, canonical=TRUE), nparam(model, canonical=TRUE))
  
  z <- NULL
  for(i in 1:nsim){
    if (verbose) {
      message(paste("#", i, " of ", nsim, ": ", sep=""),appendLF=FALSE)
    }
    
    tau <- control$SAN.tau * (if(nsim>1) (1/i-1/nsim)/(1-1/nsim) else 0)
    
    z <- ergm_SAN_slave(Clist, proposal, stats, tau, control, verbose,..., prev.run=z)

    if(z$status!=0) stop("Error in SAN.")
    
    outnw <- pending_update_network(nw,z,response=response)
    nw <-  outnw
    stats <- z$s[nrow(z$s),]
    # TODO: Subtract off a smoothing of these values?
    # Use only a bit less than the second half of the sample statistics to compute the statistic weights:
    invcov <- ginv(cov(window(z$s, start=floor(nrow(z$s)*5/8))))
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
    }else{
      if(i<nsim && isTRUE(all.equal(stats, numeric(length(stats))))){
        if(verbose) message("Target statistics matched exactly.")
        break
      }
    }
  }
  if(nsim > 1 && !only.last){
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

#' @describeIn san Sufficient statistics and other settings are
#'   inherited from the [`ergm`] fit unless overridden.
#' @export
san.ergm <- function(object, formula=object$formula, 
                     constraints=object$constraints, 
                     target.stats=object$target.stats,
                     nsim=1, basis=NULL,
                     output=c("network","edgelist","pending_update_network"),
                     only.last=TRUE,
                     control=object$control$SAN.control,
                     verbose=FALSE, ...) {
  san.formula(formula, nsim=nsim, 
              target.stats=target.stats,
              basis=basis,
              reference = object$reference,
              output=output,
              constraints=constraints,
              control=control,
              verbose=verbose, ...)
}

ergm_SAN_slave <- function(Clist,proposal,stats,tau,control,verbose,...,prev.run=NULL, nsteps=NULL, samplesize=NULL, maxedges=NULL) {
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
              as.double(c(Clist$inputs,proposal$inputs)), as.double(.deinf(tau)),
              s = as.double(c(stats, numeric(length(stats)*(samplesize-1)))),
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
              PACKAGE="ergm")
      
        # save the results (note that if prev.run is NULL, c(NULL$s,z$s)==z$s.
      z<-list(s=matrix(z$s, ncol=Clist$nstats, byrow = TRUE),
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
              as.double(c(Clist$inputs,proposal$inputs)), as.double(.deinf(tau)),
              s = as.double(c(stats, numeric(length(stats)*(samplesize-1)))),
              as.integer(samplesize),
              as.integer(nsteps), 
              newnwtails = integer(maxedges),
              newnwheads = integer(maxedges),
              newnwweights = double(maxedges),
              as.double(control$invcov),
              as.integer(verbose), 
              as.integer(maxedges),
              status = integer(1),
              PACKAGE="ergm")
      # save the results
      z<-list(s=matrix(z$s, ncol=Clist$nstats, byrow = TRUE),
              newnwtails=z$newnwtails, newnwheads=z$newnwheads, newnwweights=z$newnwweights, status=z$status, maxedges=maxedges)
    }
    
    z$s <- rbind(prev.run$s,z$s)
    
    if(z$status!=1) return(z) # Handle everything except for MCMC_TOO_MANY_EDGES elsewhere.
    
    # The following is only executed (and the loop continued) if too many edges.
    maxedges <- maxedges * 10
    if(!is.null(control$SAN.max.maxedges)){
      if(maxedges == control$SAN.max.maxedges*10) # True iff the previous maxedges exactly equaled control$SAN.max.maxedges and that was too small.
        return(z) # This will kick the too many edges problem upstairs, so to speak.
      maxedges <- min(maxedges, control$MCMC.max.maxedges)
    }
  }
}
