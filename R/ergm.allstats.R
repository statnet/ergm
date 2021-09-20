#  File R/ergm.allstats.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################
#=========================================================================
# This file contains the 2 following functions for calculating exact
# log-likelihoods
#           <ergm.allstats>
#           <ergm.exact>
#=========================================================================



##########################################################################
# The <ergm.allstats> function visits every possible network via the
# recursive algorithm in <AllStatistics.C> and tabulates the unique
# statistics.
#
# --PARAMETERS--
#   formula:  a formula of the form 'nw ~ modelterm(s)'
#   zeroobs:  whether the statistics should be shifted by the observed
#             stats (T or F); default = TRUE
#   force  :  whether to force the calculation of the statistics, despite
#             the "large" size of the network (T or F). The network is 
#             specified by 'formula' and currently "large" means an 
#             undirected network of size > 8 or directed with size > 6;
#             default = FALSE
#   maxNumChangeStatVectors: the maximum number of unique vectors of
#             statistics that may be returned; if the number of unique
#             stats vectors exceeds this count, an ERROR will occur;
#             default = 2^16
#   ...    :  extra arguments passed by <ergm.exact>; all will be ignored
#
# --RETURNED--
#   NULL:  if the network is too "large" and 'force'=FALSE, otherwise
#   a list with 2 following components:
#     weights: the proportion of networks with the statistics given in
#              the 'statmat'
#     statmat: the unique vectors of statistics, rowbound
#
##########################################################################



#' Calculate all possible vectors of statistics on a network for an ERGM
#' 
#' \code{ergm.allstats} produces a matrix of network statistics for an
#' arbitrary \code{statnet} exponential-family random graph model.  One
#' possible use for this function is to calculate the exact loglikelihood
#' function for a small network via the \code{\link{ergm.exact}} function.
#' 
#' The mechanism for doing this is a recursive algorithm, where the number of
#' levels of recursion is equal to the number of possible dyads that can be
#' changed from 0 to 1 and back again.  The algorithm starts with the network
#' passed in \code{formula}, then recursively toggles each edge twice so that
#' every possible network is visited.
#' 
#' \code{ergm.allstats} should only be used for small networks, since the
#' number of possible networks grows extremely fast with the number of nodes.
#' An error results if it is used on a directed network of more than 6 nodes or
#' an undirected network of more than 8 nodes; use \code{force=TRUE} to
#' override this error.
#' 
#' @param formula an \code{\link{formula}} object of the form \code{y ~ <model
#' terms>}, where \code{y} is a network object or a matrix that can be coerced
#' to a \code{\link[network]{network}} object.  For the details on the possible
#' \code{<model terms>}, see \code{\link{ergm-terms}}.  To create a
#' \code{\link[network]{network}} object in , use the \code{network()}
#' function, then add nodal attributes to it using the \code{\%v\%} operator if
#' necessary.
#' @param zeroobs Logical: Should the vectors be centered so that the network
#' passed in the \code{formula} has the zero vector as its statistics?
#' @param force Logical: Should the algorithm be run even if it is determined
#' that the problem may be very large, thus bypassing the warning message that
#' normally terminates the function in such cases?
#' @param maxNumChangeStatVectors Maximum possible number of distinct values of
#' the vector of statistics.  It's good to use a power of 2 for this.
#' @param \dots further arguments; not currently used.
#' @return Returns a list object with these two elements:
#' \item{weights}{integer of counts, one for each row of \code{statmat} telling
#' how many networks share the corresponding vector of statistics.}
#' \item{statmat}{matrix in which each row is a unique vector of statistics.}
#' @seealso \code{\link{ergm.exact}}
#' @keywords models
#' @examples
#' 
#' # Count by brute force all the edge statistics possible for a 7-node 
#' # undirected network
#' mynw <- network(matrix(0,7,7),dir=FALSE)
#' system.time(a <- ergm.allstats(mynw~edges))
#' 
#' # Summarize results
#' rbind(t(a$statmat),a$weights)
#' 
#' # Each value of a$weights is equal to 21-choose-k, 
#' # where k is the corresponding statistic (and 21 is 
#' # the number of dyads in an 7-node undirected network).  
#' # Here's a check of that fact:
#' as.vector(a$weights - choose(21, t(a$statmat)))
#' 
#' # Simple ergm.exact outpuf for this network.  
#' # We know that the loglikelihood for my empty 7-node network
#' # should simply be -21*log(1+exp(eta)), so we may check that
#' # the following two values agree:
#' -21*log(1+exp(.1234)) 
#' ergm.exact(.1234, mynw~edges, statmat=a$statmat, weights=a$weights)
#' 
#' @export ergm.allstats
ergm.allstats <- function(formula, zeroobs = TRUE, force = FALSE,
                          maxNumChangeStatVectors = 2^16, ...)
{
  on.exit(ergm_Cstate_clear())

  # Initialization stuff
  nw <- ergm.getnetwork(formula)
  m <- ergm_model(formula, nw, ...)

  # Check for networks that are too large.  Pretty unsophisticated check for now.
  if ((network.size(nw) > 8 && !is.directed(nw)) || (network.size(nw) > 6 && is.directed(nw))) {
    warning("Network may be too large to calculate all possible change statistics",
            immediate.=TRUE)
    if (!force) {
      cat ("Use 'force=TRUE' to proceed with this calculation.\n")
      return(NULL)
    }
  }
  
  # Calculate the statistics relative to the current network using recursive C code
  z <- .Call("AllStatistics",
             ergm_state(nw, model=m),
             # Allstats settings
             as.integer(maxNumChangeStatVectors),
          PACKAGE="ergm")
  nz <- z[[2]] > 0
  weights <- z[[2]][nz]
  statmat <- matrix(z[[1]], ncol=nparam(m,canonical=TRUE), byrow=TRUE)[nz, , drop=FALSE]

  # Add in observed statistics if they're not supposed to be zero  
  if (!zeroobs) {
    obsstats <- summary(m, nw)
    statmat <- sweep(statmat, 2, obsstats, '+')
  }

  colnames(statmat) <- param_names(m, TRUE)
  list(weights=weights, statmat=statmat)
}



##########################################################################
# The <ergm.exact> function computes the exact log-likelihood of a
# given vector of coefficients, by use of the <ergm.allstats> function.
#
# --PARAMETERS--
#   eta    :  a vector of coefficients for the model terms given in
#             'formula'
#   formula:  a formula of the form 'nw ~ modelterm(s)'
#   statmat:  the 'statmat' returned by <ergm.allstats>. This saves 
#             repeatedly calculating 'statmat' in the event that
#             <ergm.exact> is called multiple times. If using this
#             approach, 'zeroobs' should be TRUE, else the loglikelihood
#             will be shifted by the value eta %*% observed stats;
#             default=NULL
#   weights:  the 'weights' returned by <ergm.allstats>; these are used 
#             in conjuction with 'statmat' to avoid repeated calculations;
#             default=NULL
#   ...    :  additional arguments passed to <ergm.exact> that will be
#             ignored
#
# --RETURNED--
#   the exact log-likelihood at eta for the given formula IF the value
#   returned by <ergm.allstats> is non-NULL
##########################################################################



#' Calculate the exact loglikelihood for an ERGM
#' 
#' \code{ergm.exact} calculates the exact loglikelihood, evaluated at
#' \code{eta}, for the \code{statnet} exponential-family random graph model
#' represented by \code{formula}.
#' 
#' \code{ergm.exact} should only be used for small networks, since the number
#' of possible networks grows extremely fast with the number of nodes.  An
#' error results if it is used on a directed network of more than 6 nodes or an
#' undirected network of more than 8 nodes; use \code{force=TRUE} to override
#' this error.
#' 
#' In case this function is to be called repeatedly, for instance by an
#' optimization routine, it is preferable to call \code{\link{ergm.allstats}}
#' first, then pass \code{statmat} and \code{weights} explicitly to avoid
#' repeatedly calculating these objects.
#' 
#' @param eta vector of canonical parameter values at which the loglikelihood
#' should be evaluated.
#' @param formula an \code{link{formula}} object of the form \code{y ~ <model
#' terms>}, where \code{y} is a network object or a matrix that can be coerced
#' to a \code{\link[network]{network}} object.  For the details on the possible
#' \code{<model terms>}, see \code{\link{ergm-terms}}.  To create a
#' \code{\link[network]{network}} object in , use the \code{network()}
#' function, then add nodal attributes to it using the \code{\%v\%} operator if
#' necessary.
#' @param statmat if NULL, call \code{\link{ergm.allstats}} to generate all
#' possible graph statistics for the networks in this model.
#' @param weights In case \code{statmat} is not \code{NULL}, this should be the
#' vector of counts corresponding to the rows of \code{statmat}. If
#' \code{statmat} is \code{NULL}, this is generated by the call to
#' \code{\link{ergm.allstats}}.
#' @param \dots further arguments; not currently used.
#' @return Returns the value of the exact loglikelihood, evaluated at
#' \code{eta}, for the \code{statnet} exponential-family random graph model
#' represented by \code{formula}.
#' @seealso \code{\link{ergm.allstats}}
#' @keywords models
#' @examples
#' 
#' # Count by brute force all the edge statistics possible for a 7-node 
#' # undirected network
#' mynw <- network(matrix(0,7,7),dir=FALSE)
#' system.time(a <- ergm.allstats(mynw~edges))
#' 
#' # Summarize results
#' rbind(t(a$statmat),a$weights)
#' 
#' # Each value of a$weights is equal to 21-choose-k, 
#' # where k is the corresponding statistic (and 21 is 
#' # the number of dyads in an 7-node undirected network).  
#' # Here's a check of that fact:
#' as.vector(a$weights - choose(21, t(a$statmat)))
#' 
#' # Simple ergm.exact outpuf for this network.  
#' # We know that the loglikelihood for my empty 7-node network
#' # should simply be -21*log(1+exp(eta)), so we may check that
#' # the following two values agree:
#' -21*log(1+exp(.1234)) 
#' ergm.exact(.1234, mynw~edges, statmat=a$statmat, weights=a$weights)
#' 
#' @export ergm.exact
ergm.exact <- function(eta, formula, statmat=NULL, weights=NULL, ...) {
  if (is.null(statmat)) {
    alst <- ergm.allstats(formula, zeroobs=TRUE, ...)
    return(-log(alst$weights %*% exp(as.matrix(alst$statmat) %*% eta)))
  } else {
    return(-log(weights %*% exp(as.matrix(statmat) %*% eta)))
  }
}

