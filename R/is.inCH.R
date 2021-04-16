#  File R/is.inCH.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################
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

########
# The n-vector p is in the convex hull of the n-vectors in the
# Rxn matrix M iff the maximum value of z_0 + z'p equals zero for all vectors z
# satisfying z_0 + z'M \le 0.  Thus, the value returned by is.inCH depends on the
# solution of a linear program, for which the lpSolve package and its function
# lp is needed.  Letting q=(1 p')' and L = (1 M), the program is:
#
#  maximize z'q  
#  subject to z'L <= 0 
#
#  Notice that, if there exists a z that makes z'q positive, then there is
#  no maximizer since in that case Kz gives a larger value whenever K>1.  
#  For this reason, we add one additional constraint, namely, z'q <= 1.
#
#  To put this all in "standard form", we let z=a-b, where a, b, are nonnegative.
#  If we then write x=(a' b')', we obtain a new linear program:  
#
#  Minimize x'(-q' q')'
#  subject to x'(q' -q')' <= 1 and x'(L -L) <= 0 and x >= 0
#  ...and if the minimum is strictly negative, return FALSE because the point
#  is not in the CH in that case.

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
#' Letting \eqn{q=(1 p')'} and \eqn{L = (1 M)}, if the maximum value of
#' \eqn{z'q} for all \eqn{z} such that \eqn{z'L \le 0} equals zero (the maximum
#' must be at least zero since z=0 gives zero), then there is no separating
#' hyperplane and so \code{p} is contained in the convex hull of the rows of
#' \code{M}.  So the question of interest becomes a constrained optimization
#' problem.
#' 
#' Solving this problem relies on the package \code{lpSolve} to solve a linear
#' program.  We may put the program in "standard form" by writing \eqn{z=a-b},
#' where \eqn{a} and \eqn{b} are nonnegative vectors.  If we write \eqn{x=(a'
#' b')'}, we obtain the linear program given by:
#' 
#' Minimize \eqn{(-q' q')x} subject to \eqn{x'(L -L) \le 0} and \eqn{x \ge 0}.
#' One additional constraint arises because whenever any strictly negative
#' value of \eqn{(-q' q')x} may be achieved, doubling \eqn{x} arbitrarily many
#' times makes this value arbitrarily large in the negative direction, so no
#' minimizer exists.  Therefore, we add the constraint \eqn{(q' -q')x \le 1}.
#' 
#' This function is used in the "stepping" algorithm of Hummel et al (2012).
#' 
#' @param p A \eqn{d}-dimensional vector or a matrix with \eqn{d} columns
#' @param M An \eqn{r} by \eqn{d} matrix.  Each row of \code{M} is a
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
#' @export is.inCH
is.inCH <- function(p, M, verbose=FALSE, ...) { # Pass extra arguments directly to LP solver

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
  L <- cbind(1, M)
  lprec <- make.lp(n, d+1)
  for(k in seq_len(d+1)) set.column(lprec, k, L[,k])
  set.constr.type(lprec, rep.int(2L, n)) # 2 = ">="
  set.rhs(lprec,  numeric(n))
  set.bounds(lprec, lower = rep.int(-1, d+1L), upper = rep.int(1, d+1L))

  for(i in seq_len(nrow(p))){ # Iterate over test points.

    # Set the objective function in terms of p and solve the problem.
    set.objfn(lprec, c(1, p[i,]))
    solve(lprec)

    # If the objective function (min) is not zero, the point p[i,] is not in the CH of M.
    if(get.objective(lprec)){
      if(verbose>1) message(sprintf("is.inCH: iter= %d, outside hull.",i))
      return(FALSE)
    }
  }

  if(verbose>1) message(sprintf("is.inCH: iter= %d, inside hull.",i))
  return(TRUE) # If all points passed the test, return TRUE.
}
