#  File R/ergm.san.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2018 Statnet Commons
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
#' @param sequential Logical: If TRUE, the returned draws always use the prior
#' draw as the starting network; if FALSE, they always use the original
#' network.
#' @param control A list of control parameters for algorithm tuning; see
#' \code{\link{control.san}}.
#' @param verbose Logical: If TRUE, print out more detailed information as the
#' simulation runs.
#' @param \dots Further arguments passed to other functions.
#' @export
san.formula <- function(object, response=NULL, reference=~Bernoulli, constraints=~., target.stats=NULL,
                        nsim=1, basis=NULL,
                        sequential=TRUE,
                        control=control.san(),
                        verbose=FALSE, ...) {
  check.control.class("san", "san")
  control.toplevel(...,myname="san")

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
  nw <- as.network(nw)
  if(is.null(target.stats)){
    stop("You need to specify target statistic via",
         " the 'target.stats' argument")
  }
  if(!is.network(nw)){
    stop("A network object on the LHS of the formula ",
         "must be given")
  }

# model <- ergm_model(formula, nw, drop=control$drop)
  model <- ergm_model(formula, nw, response=response, term.options=control$term.options)
  Clist <- ergm.Cprepare(nw, model, response=response)
  fd <- ergm.design(nw, verbose=verbose)
  
  verb <- match(verbose,
                c("FALSE","TRUE", "very"), nomatch=1)-1
  proposal<-ergm_proposal(constraints,arguments=control$SAN.prop.args,nw=nw,weights=control$SAN.prop.weights, class="c",reference=reference,response=response)
# if(is.null(control$coef)) {
#   warning("No parameter values given, using the MPLE for the passed network.\n\t")
# }
# control$coef <- c(control$coef[1],rep(0,Clist$nstats-1))
  
  if (verb) {
    message(paste("Starting ",nsim," MCMC iteration", ifelse(nsim>1,"s",""),
        " of ", control$SAN.burnin+control$SAN.interval*(nsim-1), 
        " steps", ifelse(nsim>1, " each", ""), ".", sep=""))
  }

  for(i in 1:nsim){
    Clist <- ergm.Cprepare(nw, model,response=response)
#   fd <- ergm.design(nw, verbose=verbose)
    maxedges <- max(control$SAN.init.maxedges, Clist$nedges)
    if (verb) {
       message(paste("#", i, " of ", nsim, ": ", sep=""),appendLF=FALSE)
     }

    if(is.null(control$coef)) {
      if(reference==~Bernoulli){
        fit <- suppressWarnings(suppressMessages(try(ergm.mple(nw=nw, fd=fd, 
                         control=control, proposal=proposal,
                         m=model, verbose=verbose, ...))))
        control$coef <- if(inherits(fit, "try-error")) rep(0,nparam(model,canonical=TRUE)) else fit$coef
        if(is.null(control$invcov)) { control$invcov <- fit$covar }
      }else{
        control$coef<-rep(0,nparam(model,canonical=TRUE))
        if(is.null(control$invcov)) control$invcov <- diag(length(control$coef))
      }
    }else{
      if(is.null(control$invcov)) { control$invcov <- diag(length(control$coef)) }
    }
    eta0 <- ifelse(is.na(control$coef), 0, control$coef)
    
    netsumm<-summary(model,nw,response=response)
    target.stats <- vector.namesmatch(target.stats, names(netsumm))
    
    stats <- matrix(netsumm-target.stats,
                    ncol=Clist$nstats,byrow=TRUE,nrow=nsim)
    tau <- rep(control$SAN.tau,length=length(eta0))
#
#   Check for truncation of the returned edge list
#
    repeat{
      nedges <- c(Clist$nedges,0)
      tails <- Clist$tails
      heads <- Clist$heads
      weights <- Clist$weights
      # *** don't forget to pass in tails before heads now.

      if(is.null(Clist$weights)){
        z <- .C("SAN_wrapper",
                as.integer(length(nedges)), as.integer(nedges),
                as.integer(tails), as.integer(heads),
                as.integer(Clist$n),
                as.integer(Clist$dir), as.integer(Clist$bipartite),
                as.integer(Clist$nterms), 
                as.character(Clist$fnamestring),
                as.character(Clist$snamestring), 
                as.character(proposal$name),
                as.character(proposal$pkgname),
                as.double(c(Clist$inputs,proposal$inputs)),
                as.double(.deinf(eta0)),
                as.double(.deinf(tau)),
                as.integer(1), # "samplesize"
                s = as.double(stats),
                as.integer(if(i==1 | !sequential) control$SAN.burnin else control$SAN.interval), as.integer(control$SAN.interval), 
                newnwtails = integer(maxedges),
                newnwheads = integer(maxedges), 
                as.double(control$invcov),
                as.integer(verb),
                as.integer(proposal$arguments$constraints$bd$attribs), 
                as.integer(proposal$arguments$constraints$bd$maxout), as.integer(proposal$arguments$constraints$bd$maxin),
                as.integer(proposal$arguments$constraints$bd$minout), as.integer(proposal$arguments$constraints$bd$minin),
                as.integer(proposal$arguments$constraints$bd$condAllDegExact),
                as.integer(length(proposal$arguments$constraints$bd$attribs)), 
                as.integer(maxedges),
                status = integer(1),
                PACKAGE="ergm")
      }else{
        z <- .C("WtSAN_wrapper",
                as.integer(length(nedges)), as.integer(nedges),
                as.integer(tails), as.integer(heads), as.double(weights),
                as.integer(Clist$n),
                as.integer(Clist$dir), as.integer(Clist$bipartite),
                as.integer(Clist$nterms), 
                as.character(Clist$fnamestring),
                as.character(Clist$snamestring), 
                as.character(proposal$name),
                as.character(proposal$pkgname),
                as.double(c(Clist$inputs,proposal$inputs)),
                as.double(.deinf(eta0)),
                as.double(.deinf(tau)),
                as.integer(1), # "samplesize"
                s = as.double(stats),
                as.integer(if(i==1 | !sequential) control$SAN.burnin else control$SAN.interval), as.integer(control$SAN.interval), 
                newnwtails = integer(maxedges),
                newnwheads = integer(maxedges),
                newnwweights = double(maxedges), 
                as.double(control$invcov),
                as.integer(verb),
                as.integer(proposal$arguments$constraints$bd$attribs), 
                as.integer(proposal$arguments$constraints$bd$maxout), as.integer(proposal$arguments$constraints$bd$maxin),
                as.integer(proposal$arguments$constraints$bd$minout), as.integer(proposal$arguments$constraints$bd$minin),
                as.integer(proposal$arguments$constraints$bd$condAllDegExact),
                as.integer(length(proposal$arguments$constraints$bd$attribs)), 
                as.integer(maxedges), 
                status = integer(1),
                PACKAGE="ergm")
      }

      if(z$status==1) maxedges <- 5*maxedges
      else if(z$status==0) break
      else stop("Error in SAN C code.")
    }
    #
    #   Next update the network to be the final (possibly conditionally)
    #   simulated one
    #
    out.list[[i]] <- newnw.extract(nw, z, output=control$network.output, response=response)
    out.mat <- rbind(out.mat,z$s[(Clist$nstats+1):(2*Clist$nstats)])
    if(sequential){
      nw <-  as.network.uncompressed(out.list[[i]])
    }
  }
  if(nsim > 1){
    out.list <- list(formula = formula, networks = out.list, 
                     stats = out.mat, coef=control$coef)
    class(out.list) <- "network.list"
  }else{
    out.list <- out.list[[1]]
  }
  return(out.list)
}

#' @describeIn san Sufficient statistics and other settings are
#'   inherited from the [`ergm`] fit unless overridden.
#' @export
san.ergm <- function(object, formula=object$formula, 
                     constraints=object$constraints, 
                     target.stats=object$target.stats,
                     nsim=1, basis=NULL,
                     sequential=TRUE, 
                     control=object$control$SAN.control,
                     verbose=FALSE, ...) {
  if(is.null(control$coef)) control$coef <- coef(object)
  san.formula(formula, nsim=nsim, 
              target.stats=target.stats,
              basis=basis,
              reference = object$reference,
              sequential=sequential,
              constraints=constraints,
              control=control,
              verbose=verbose, ...)
}

