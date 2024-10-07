#  File R/is.inCH.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################

warning_once <- once(warning)

#' Identify the position of a point relative to the  convex hull of a set of points
#' 
#' This function uses linear programming to find the value by which
#' vector `p` needs to be scaled towards or away from vector `m` in
#' order for `p` to be on the boundary of the convex hull of rows of
#' `M`. If `p` is a matrix, a value that scales all rows of `p` into
#' the convex hull of `M` is found.
#'
#' @note This is a successor to the deprecated function `is.inCH()`,
#'   which was originally written for the "stepping" algorithm of
#'   \insertCite{HuHu12i;textual}{ergm}. See the updated of
#'   \insertCite{KrKu23l;textual}{ergm} for detailed discussion of algorithms
#'   used in `is.inCH()` and `shrink_into_CH()`.
#'
#' @param p a \eqn{d}-dimensional vector or a matrix with \eqn{d}
#'   columns.
#' @param M an \eqn{n} by \eqn{d} matrix.  Each row of \code{M} is a
#'   \eqn{d}-dimensional vector.
#' @param m a \eqn{d}-dimensional vector specifying the value towards
#'   which to shrink; must be in the interior of the convex hull of
#'   \eqn{M}, and defaults to its centroid (column means).
#' @template verbose
#' @param max_run if there are no decreases in step length in this
#'   many consecutive test points, conclude that diminishing returns
#'   have been reached and finish.
#' @param \dots arguments passed directly to linear program solver.
#' @param solver a character string selecting which solver to use; by
#'   default, tries `Rglpk`'s but falls back to `lpSolveAPI`'s.
#'
#' @return The scaling factor described above is
#'   returned. `shrink_into_CH() >= 1` indicates that all points in
#'   `p` are in the convex hull of `M`.
#'
#' @references \insertAllCited{}
#'
#' \url{https://www.cs.mcgill.ca/~fukuda/soft/polyfaq/node22.html}
#'
#' @keywords internal
#' @export
shrink_into_CH <- function(p, M, m = NULL, verbose = FALSE, max_run = nrow(M), ..., solver = c("glpk", "lpsolve")) { # Pass extra arguments directly to LP solver
  solver <- match.arg(solver)
  verbose <- max(0, min(verbose, 4))

  if(solver == "glpk" && !requireNamespace("Rglpk", quietly=TRUE)){
    warning_once(sQuote("glpk"), " selected as the solver, but package ", sQuote("Rglpk"), " is not available; falling back to ", sQuote("lpSolveAPI"), ". This should be fine unless the sample size and/or the number of parameters is very big.", immediate.=TRUE, call.=FALSE)
    solver <- "lpsolve"
  }

  if(is.null(dim(p))) p <- rbind(p)

  if (!is.matrix(M))
    stop("Second argument must be a matrix.")
  if ((d <- ncol(p)) != ncol(M))
    stop("Number of columns in matrix (2nd argument) is not equal to dimension ",
         "of first argument.")

  NVL(m) <- colMeans(M)
  p <- sweep_cols.matrix(p, m)
  np <- nrow(p)
  M <- sweep_cols.matrix(M, m)

  if((n <- nrow(M)) == 1L){
    for(i in seq_len(np)){
      if(!isTRUE(all.equal(p[i,], M, check.attributes = FALSE))) return(0)
    }
    return(1)
  }

  # Minimise: p'z
  # Constrain: Mz >= -1. No further constraints!
  dir <- rep.int(">=", n)
  rhs <- rep.int(-1, n)
  lb <- rep.int(-Inf, d)

  if(solver == "lpsolve"){
    #' @importFrom lpSolveAPI make.lp set.column set.objfn set.constr.type set.rhs set.bounds get.objective
    setup.lp <- function(){
      # Set up the optimisation problem: the following are common for all rows of p.
      lprec <- make.lp(n, d)
      for(k in seq_len(d)) set.column(lprec, k, M[, k])
      set.constr.type(lprec, dir)
      set.rhs(lprec,  rhs)
      # By default, z are bounded >= 0. We need to remove these bounds.
      set.bounds(lprec, lower=lb)
      lp.control(lprec, verbose=c("important","important","important","normal","detailed")[min(max(verbose+1,0),5)], ...)
      lprec
    }
    lprec <- setup.lp()
  }else{
    # Rglpk prefers this format.
    M <- slam::as.simple_triplet_matrix(M)
  }

  if (verbose >= 2) message("Iterating over at most ", np, " test points:")
  g <- Inf
  run <- 0L
  for (i in seq_len(np)) { # Iterate over test points.
    message(i, " ", appendLF=FALSE)
    if (all(abs((x <- p[i,])) <= sqrt(.Machine$double.eps))) next # Test point is at centroid. TODO: Allow the user to specify tolerance?

    if(solver == "lpsolve"){
      # Keep trying until results are satisfactory.
      #
      # flag meanings:
      # -1      : dummy value, just starting out
      #  0 or 11: Good (either 0 or some negative value)
      #  1 or  7: Timeout
      #  3      : Unbounded: sometimes happens and solved by reinitializing
      #   others: probably nothing good, but don't know how to handle
      flag <- -1
      while(flag%in%c(-1,1,7,3)){
        set.objfn(lprec, x)
        flag <- solve(lprec)
        if(flag %in% c(1,7)){ # Timeout
          timeout <- timeout * 2 # Increase timeout, in case it's just a big problem.
          z <- rnorm(1) # Shift target and test set by the same constant.
          p <- p + z
          M <- M + z
          lprec <- setup.lp() # Reinitialize
        }else if(flag == 3){ # Unbounded
          lprec <- setup.lp() # Reinitialize
        }
      }
      o <- get.objective(lprec)
    }else{
      o <- Rglpk::Rglpk_solve_LP(x, M, dir, rhs, list(lower=list(ind=seq_len(d), val=lb)), control=list(..., verbose=max(0,verbose-3)))$optimum
    }

    g.prev <- g
    g <- min(g, abs(-1/o)) # abs() guards against optimum being numerically equivalent to 0 with -1/0 = -Inf.

    if (g < g.prev){
      if (verbose >= 3) message("|", sprintf("%0.4f", g), "| ", appendLF = FALSE)
      run <- 0L
    } else {
      run <- run + 1L
      if (run >= max_run) break
    }
  }
  if(verbose >= 3) message("")

  g
}
