#  File ergm/R/ergm.allstats.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
##########################################################################
# The <ergm.allstats> function visits every possible network via the
# recursive algorithm in <AllStatistics.C> and tabulates the unique
# statistics.
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
##########################################################################

ergm.exact <- function(eta, formula, statmat=NULL, weights=NULL, ...) {
  if (is.null(statmat)) {
    alst <- ergm.allstats(formula, zeroobs=TRUE, ...)
    return(-log(alst$weights %*% exp(as.matrix(alst$statmat) %*% eta)))
  } else {
    return(-log(weights %*% exp(as.matrix(statmat) %*% eta)))
  }
}

