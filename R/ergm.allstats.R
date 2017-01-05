#  File R/ergm.allstats.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
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

ergm.allstats <- function(formula, zeroobs = TRUE, force = FALSE,
                          maxNumChangeStatVectors = 2^16, ...)
{
  # Initialization stuff
  nw <- ergm.getnetwork(formula)
  model <- ergm.getmodel(formula, nw, initialfit=TRUE)
  Clist <- ergm.Cprepare(nw, model)

  # Check for networks that are too large.  Pretty unsophisticated check for now.
  if ((Clist$n > 8 && !Clist$dir) || (Clist$n > 6 && Clist$dir)) {
    warning("Network may be too large to calculate all possible change statistics",
            immediate.=TRUE)
    if (!force) {
      cat ("Use 'force=TRUE' to proceed with this calculation.\n")
      return(NULL)
    }
  }
  
  # Calculate the statistics relative to the current network using recursive C code
  z <- .C("AllStatistics",
          as.integer(Clist$tails),
          as.integer(Clist$heads),
          as.integer(Clist$nedges),
          as.integer(Clist$n),
          as.integer(Clist$dir),
          as.integer(Clist$bipartite),
          as.integer(Clist$nterms),
          as.character(Clist$fnamestring),
          as.character(Clist$snamestring),
          as.double(Clist$inputs),
          statmat = double(Clist$nstats * maxNumChangeStatVectors),
          weights = integer(maxNumChangeStatVectors),
          as.integer(maxNumChangeStatVectors),
          PACKAGE="ergm")
  nz <- z$weights > 0
  weights <- z$weights[nz]
  statmat <- matrix(z$statmat, ncol=Clist$nstats, byrow=TRUE)[nz, , drop=FALSE]

  # Add in observed statistics if they're not supposed to be zero  
  if (!zeroobs) {
    obsstats <- summary(formula)
    statmat <- sweep(statmat, 2, obsstats, '+')
  }
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

ergm.exact <- function(eta, formula, statmat=NULL, weights=NULL, ...) {
  if (is.null(statmat)) {
    alst <- ergm.allstats(formula, zeroobs=TRUE, ...)
    return(-log(alst$weights %*% exp(as.matrix(alst$statmat) %*% eta)))
  } else {
    return(-log(weights %*% exp(as.matrix(statmat) %*% eta)))
  }
}

