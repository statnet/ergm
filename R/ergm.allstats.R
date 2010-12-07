# ergm.allstats produces matrix of network statistics for an arbitrary 
# statnet model by an algorithm that starts with the network passed in the
# formula, then recursively toggling each edge two times so that every 
# possible network is visited.  
#
# This should only be attempted for smallish networks, since the number of
# possible networks grows extremely fast with the number of nodes.
#
# output of this function:  A matrix of statistics in which each row is
# one vector of statistics, and a vector of integer counts, one for each row
# telling how many networks share those particular statistics.
ergm.allstats <- function(formula, zeroobs = TRUE, force = FALSE,
                          maxNumChangeStatVectors = 2^16, ...)
{
  # Initialization stuff
  nw <- ergm.getnetwork(formula)
  model <- ergm.getmodel(formula, nw, drop=FALSE, initialfit=TRUE)
  Clist <- ergm.Cprepare(nw, model)
  Clist.miss <- ergm.design(nw, model, verbose=verbose)

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
          as.integer(Clist$heads),
          as.integer(Clist$tails),
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

# Use ergm.allstats output to calculate the exact loglikelihood at eta
# for a particular formula.  If repeated calls of this function will be
# needed, it's best to call ergm.allstats once first, then pass in the
# statmatrix and corresponding weights explicitly.  This should be done
# with zeroobs = TRUE, or else the loglikelihood will be shifted by
# the value eta %*% observed stats.
ergm.exact <- function(eta, formula, statmat=NULL, weights=NULL, ...) {
  if (is.null(statmat)) {
    alst <- ergm.allstats(formula, zeroobs=TRUE, ...)
    return(-log(alst$weights %*% exp(as.matrix(alst$statmat) %*% eta)))
  } else {
    return(-log(weights %*% exp(as.matrix(statmat) %*% eta)))
  }
}

