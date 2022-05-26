#  File R/is.inCH.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
###############################################################################
# The <is.inCH> function determines whether a vector p is in the convex hull
# of the vectors M
#
# --PARAMETERS--
#   p:  a vector of length n
#   M:  a q by n matrix 
#
# --RETURNED--
#   x: TRUE if p is in the CH of the points M
#      FALSE othewise
#
###############################################################################

#' Determine whether a vector is in the closure of the convex hull of some
#' sample of vectors
#' 
#' \code{is.inCH()} returns \code{TRUE} if and only if \code{p} is contained in
#' the convex hull of the points given as the rows of \code{M}. If \code{p} is
#' a matrix, each row is tested individually, and \code{TRUE} is returned if
#' all rows are in the convex hull.
#'
#' `is.inCH()` was originally written for the "stepping" algorithm of
#' Hummel et al (2012). See Krivitsky, Kuvelkar, and Hunter (2022) for
#' detailed discussion of algorithms used in `is.inCH()` and
#' `shrink_into_CH()`.
#'
#' @note [is.inCH()] has been deprecated in favour of
#'   [shrink_into_CH()], which returns the optimal step length instead
#'   of a yes-or-no test. In general, `shrink_into_CH(...)>=1` is
#'   equivalent to `is.inCH(...).
#' 
#' @param p A \eqn{d}-dimensional vector or a matrix with \eqn{d} columns
#' @param M An \eqn{n} by \eqn{d} matrix.  Each row of \code{M} is a
#' \eqn{d}-dimensional vector.
#' @template verbose
#' @param \dots arguments passed directly to linear program solver
#' @return Logical, telling whether \code{p} is (or all rows of \code{p} are)
#' in the closed convex hull of the points in \code{M}.
#' @references \itemize{ \item
#' \url{https://www.cs.mcgill.ca/~fukuda/soft/polyfaq/node22.html}
#' 
#' \item Hummel, R. M., Hunter, D. R., and Handcock, M. S. (2012), Improving
#' Simulation-Based Algorithms for Fitting ERGMs, *Journal of Computational and
#' Graphical Statistics*, 21: 920-939.
#'
#' \item Krivitsky, P. N., Kuvelkar, A. R., and Hunter,
#' D. R. (2022). Likelihood-based Inference for Exponential-Family
#' Random Graph Models via Linear Programming. *arXiv preprint*
#' arXiv:2202.03572. \url{https://arxiv.org/abs/2202.03572}
#'
#' }
#' @keywords internal
#' @export is.inCH
is.inCH <- function(p, M, verbose=FALSE, ...) { # Pass extra arguments directly to LP solver
  .Deprecate_once("shrink_into_CH()")

  verbose <- max(0, min(verbose, 4))

  if(is.null(dim(p))) p <- rbind(p)

  if (!is.matrix(M)) 
    stop("Second argument must be a matrix.")
  if ((d <- ncol(p)) != ncol(M))
    stop("Number of columns in matrix (2nd argument) is not equal to dimension ",
         "of first argument.")

  if((n <- nrow(M)) == 1L){
    for(i in seq_len(nrow(p))){
      if(!isTRUE(all.equal(p[i,], M, check.attributes = FALSE))) return(FALSE)
    }
    return(TRUE)
  }

  #' @importFrom lpSolveAPI make.lp set.column set.objfn set.constr.type set.rhs set.bounds get.objective

  # Set up the optimisation problem: the following are common for all rows of p.

  timeout <- 1

  setup.lp <- function(){
    L <- cbind(1, M)
    lprec <- make.lp(n, d+1)
    for(k in seq_len(d+1)) set.column(lprec, k, L[,k])
    set.constr.type(lprec, rep.int(2L, n)) # 2 = ">="
    set.rhs(lprec,  numeric(n))
    set.bounds(lprec, lower = rep.int(-1, d+1L), upper = rep.int(1, d+1L))
    lp.control(lprec, break.at.value = -.Machine$double.eps, verbose=c("important","important","important","normal","detailed")[min(max(verbose+1,0),5)], timeout=timeout)
    lprec
  }
  lprec <- setup.lp()

  for(i in seq_len(nrow(p))){ # Iterate over test points.

    # Keep trying until results are satisfactory.
    #
    # flag meanings:
    # -1      : dummy value, just starting out
    #  0 or 11: Good (either 0 or some negative value)
    #  1 or  7: Timeout
    #   others: probably nothing good, but don't know how to handle
    flag <- -1
    while(flag%in%c(-1,1,7)){
      # Set the objective function in terms of p and solve the problem.
      set.objfn(lprec, c(1, p[i,]))
      flag <- solve(lprec)
      if(flag%in%c(1,7)){ # Timeout
        timeout <- timeout * 2 # Increase timeout, in case it's just a big problem.
        z <- rnorm(1) # Shift target and test set by the same constant.
        p <- p + z
        M <- M + z
        lprec <- setup.lp() # Reinitialize
      }
    }

    # If the objective function (min) is not zero, the point p[i,] is not in the CH of M.
    if(get.objective(lprec) < 0){
      if(verbose > 1) message(sprintf("is.inCH: iter = %d, outside hull.",i))
      return(FALSE)
    }else if(verbose > 2 && nrow(p) > 1) message(sprintf("is.inCH: iter = %d, inside hull.", i))

  }

  if(verbose > 1) message("is.inCH: all test points inside hull.")
  return(TRUE) # If all points passed the test, return TRUE.
}


warning_once <- once(warning)

#' @rdname is.inCH
#' @description `shrink_into_CH()` returns the coefficient by which rows of `p` can be scaled towards or away from point `m` in order for all of them to be in the convex hull of `M` or on its boundary.
#' @param solver A character string selecting which solver to use; by default, tries `Rglpk`'s but falls back to `lpSolveAPI`'s.
#' @export
shrink_into_CH <- function(p, M, m = NULL, verbose=FALSE, ..., solver = c("glpk", "lpsolve")) { # Pass extra arguments directly to LP solver
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
    # Set up the optimisation problem: the following are common for all rows of p.
    lprec <- make.lp(n, d)
    for(k in seq_len(d)) set.column(lprec, k, M[, k])
    set.constr.type(lprec, dir)
    set.rhs(lprec,  rhs)
    # By default, z are bounded >= 0. We need to remove these bounds.
    set.bounds(lprec, lower=lb)
    lp.control(lprec, verbose=c("important","important","important","normal","detailed")[min(max(verbose+1,0),5)], ...)
  }else{
    # Rglpk prefers this format.
    M <- slam::as.simple_triplet_matrix(M)
  }

  if (verbose >= 2) message("Iterating over ", np, " test points.")
  g <- Inf
  for (i in seq_len(np)) { # Iterate over test points.
    if (verbose >= 3) message("Test point ", i)
    if (all(abs((x <- p[i,])) <= sqrt(.Machine$double.eps))) next # Test point is at centroid. TODO: Allow the user to specify tolerance?

    if(solver == "lpsolve"){
      set.objfn(lprec, x)
      solve(lprec)
      o <- get.objective(lprec)
    }else{
      o <- Rglpk::Rglpk_solve_LP(x, M, dir, rhs, list(lower=list(ind=seq_len(d), val=lb)), control=list(..., verbose=max(0,verbose-3)))$optimum
    }

    g <- min(g, abs(-1/o)) # abs() guards against optimum being numerically equivalent to 0 with -1/0 = -Inf.

    if (verbose >= 3) message("Step length is now ", g, ".")
  }

  g
}
