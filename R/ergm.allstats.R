#  File R/ergm.allstats.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################
#=========================================================================
# This file contains the 2 following functions for calculating exact
# log-likelihoods
#           <ergm.allstats>
#           <ergm.exact>
#=========================================================================


#' Calculate all possible vectors of statistics on a network for an ERGM
#' 
#' \code{ergm.allstats} calculates the sufficient statistics of an
#' ERGM over the network's sample space.
#' 
#' The mechanism for doing this is a recursive algorithm, where the number of
#' levels of recursion is equal to the number of possible dyads that can be
#' changed from 0 to 1 and back again.  The algorithm starts with the network
#' passed in \code{formula}, then recursively toggles each edge twice so that
#' every possible network is visited.
#' 
#' `ergm.allstats()` and `ergm.exact()` should only be used for small
#' networks, since the number of possible networks grows extremely
#' fast with the number of nodes.  An error results if it is used on a
#' network with more than 31 free dyads, which corresponds to a
#' directed network of more than 6 nodes or an undirected network of
#' more than 8 nodes; use \code{force=TRUE} to override this error.
#'
#' @param formula,constraints An ERGM formula and
#'   (optionally) a constraint specification formulas. See
#'   [ergm()]. This function supports only dyad-independent
#'   constraints.
#' @param zeroobs Logical: Should the vectors be centered so that the network
#' passed in the \code{formula} has the zero vector as its statistics?
#' @param force Logical: Should the algorithm be run even if it is determined
#' that the problem may be very large, thus bypassing the warning message that
#' normally terminates the function in such cases?
#' @param \dots further arguments, passed to [ergm_model()].
#' @return `ergm.allstats()` returns a list object with these two elements:
#' \item{weights}{integer of counts, one for each row of \code{statmat} telling
#' how many networks share the corresponding vector of statistics.}
#' \item{statmat}{matrix in which each row is a unique vector of statistics.}
#' @keywords models
#' @examples
#' 
#' # Count by brute force all the edge statistics possible for a 7-node 
#' # undirected network
#' mynw <- network.initialize(7, dir = FALSE)
#' system.time(a <- ergm.allstats(mynw~edges))
#' 
#' # Summarize results
#' rbind(t(a$statmat), .freq. = a$weights)
#' 
#' # Each value of a$weights is equal to 21-choose-k, 
#' # where k is the corresponding statistic (and 21 is 
#' # the number of dyads in an 7-node undirected network).  
#' # Here's a check of that fact:
#' as.vector(a$weights - choose(21, t(a$statmat)))
#'
#' # Dyad-independent constraints are also supported:
#' system.time(a <- ergm.allstats(mynw~edges, constraints = ~fixallbut(cbind(1:2,2:3))))
#' rbind(t(a$statmat), .freq. = a$weights)
#' 
#' @export ergm.allstats
ergm.allstats <- function(formula, constraints=~., zeroobs = TRUE, force = FALSE, ...)
{
  on.exit(ergm_Cstate_clear())
  on.exit(allstats_workspace_clear(), add=TRUE)

  # Initialization stuff
  nw <- ergm.getnetwork(formula)
  tmp <- .handle.auto.constraints(nw, constraints)
  nw <- tmp$nw; conterms <- tmp$conterms
  conlist <- ergm_conlist(conterms, nw)
  if(!is.dyad.independent(conlist)) stop("Only dyadic constraints are supported.")
  fd <- as.rlebdm(conlist)
  m <- ergm_model(formula, nw, ...)

  # Check for networks that are too large (based on the number of free dyads).
  if (sum(fd) > 31) { # TODO: Make tunable.
    if(force) warning("Network may be too large to calculate all possible change statistics.",
                      immediate.=TRUE)
    else stop("Network may be too large to calculate all possible change statistics. Use ", sQuote("force=TRUE")," to proceed with this calculation.")
  }
  
  # Calculate the statistics relative to the current network using recursive C code
  z <- .Call("AllStatistics",
             ergm_state(nw, model=m),
             # Allstats settings
             as.double(to_ergm_Cdouble(fd)),
             PACKAGE="ergm")
  weights <- z$weights
  statmat <- z$stats

  # Add in observed statistics if they're not supposed to be zero  
  if (!zeroobs) {
    obsstats <- summary(m, nw)
    statmat <- statmat + obsstats
  }
  statmat <- t(statmat)

  colnames(statmat) <- param_names(m, TRUE)
  list(weights=weights, statmat=statmat)
}

allstats_workspace_clear <- function(){
  .Call("allstats_workspace_free", PACKAGE="ergm")
}

#' @rdname ergm.allstats
#'
#' @description \code{ergm.exact()} uses \code{ergm.allstats()} to calculate the exact loglikelihood, evaluated at
#' \code{eta}.
#'
#' @details
#' In case \code{ergm.exact()} is to be called repeatedly, for instance by an
#' optimization routine, it is preferable to call `ergm.allstats()`
#' first, then pass \code{statmat} and \code{weights} explicitly to avoid
#' repeatedly calculating these objects.
#' 
#' @param eta vector of canonical parameter values at which the loglikelihood
#' should be evaluated.
#' @param statmat,weights outputs from `ergm.allstats()`: if passed, used in lieu of running it.
#' @return `ergm.exact()` returns the exact value of the loglikelihood, evaluated at
#' \code{eta}.
#' @keywords models
#' @examples
#' 
#' # Simple ergm.exact output for this network.
#' # We know that the loglikelihood for my empty 7-node network
#' # should simply be -21*log(1+exp(eta)), so we may check that
#' # the following two values agree:
#' -21*log(1+exp(.1234)) 
#' ergm.exact(.1234, mynw~edges, statmat=a$statmat, weights=a$weights)
#' 
#' @export ergm.exact
ergm.exact <- function(eta, formula, constraints = ~., statmat=NULL, weights=NULL, ...) {
  if (is.null(statmat)) {
    alst <- ergm.allstats(formula, constraints, zeroobs=TRUE, ...)
    statmat <- alst$statmat; weights <- alst$weights
  }

  -log(weights %*% exp(as.matrix(statmat) %*% eta))
}

