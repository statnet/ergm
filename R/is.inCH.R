#  File R/is.inCH.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
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
#' \code{is.inCH} returns \code{TRUE} if and only if \code{p} is contained in
#' the convex hull of the points given as the rows of \code{M}. If \code{p} is
#' a matrix, each row is tested individually, and \code{TRUE} is returned if
#' all rows are in the convex hull.
#' 
#' The \eqn{d}-vector \code{p} is in the convex hull of the \eqn{d}-vectors
#' forming the rows of \code{M} if and only if there exists no separating
#' hyperplane between \code{p} and the rows of \code{M}.  This condition may be
#' reworded as follows:
#' 
#' Letting \eqn{q=(1 p')'} and \eqn{L = (1 M)}, if the minimum value
#' of \eqn{z'q} for all \eqn{z} such that \eqn{z'L \ge 0} equals zero
#' (the minimum must be at most zero since z=0 gives zero), then there
#' is no separating hyperplane and so \code{p} is contained in the
#' convex hull of the rows of \code{M}. So the question of interest
#' becomes a constrained optimization problem.
#'
#' Lastly, in the event of such a hyperplane existing, one can make
#' the objective function arbitrarily low by multiplying \eqn{z} by a
#' large positive constant. To prevent it from running away, we
#' constrain the elements of \eqn{z} to be between \eqn{-1} and
#' \eqn{+1}.
#' 
#' Solving this problem relies on the package \pkg{lpSolveAPI} to solve a linear
#' program.
#' 
#' This function is used in the "stepping" algorithm of Hummel et al (2012).
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
#' Simulation-Based Algorithms for Fitting ERGMs, Journal of Computational and
#' Graphical Statistics, 21: 920-939. }
#'
#' @keywords internal
#' @export is.inCH
is.inCH <- function(p, M, verbose=FALSE, ...) { # Pass extra arguments directly to LP solver
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
