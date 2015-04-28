%  File man/is.inCH.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2015 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{is.inCH}
\alias{is.inCH}
\title{Determine whether a vector is in the closure of the convex hull of 
some sample of vectors}
\description{
  \code{is.inCH} returns \code{TRUE} if and only if 
  \code{p} is contained in the convex hull of the points given as
  the rows of \code{M}.
}
\usage{
is.inCH(p, M, \dots)
}
\arguments{
  \item{p}{A \eqn{d}-dimensional vector}
  \item{M}{An \eqn{r} by \eqn{d} matrix.  Each row of \code{M} is a 
  \eqn{d}-dimensional vector.}
  \item{\dots}{arguments passed directly to linear program solver}
}
\details{
  The \eqn{d}-vector \code{p} is in the convex hull of the \eqn{d}-vectors
  forming the rows of \code{M} if and only if there exists no separating 
  hyperplane between \code{p} and the rows of \code{M}.  This condition
  may be reworded as follows: 
  
  Letting \eqn{q=(1 p')'} and \eqn{L = (1 M)}, if the maximum value of
  \eqn{z'q} for all \eqn{z} such that \eqn{z'L \le 0} equals zero
  (the maximum must be at least zero since z=0 gives zero), then there is no 
  separating hyperplane and so \code{p} is contained in the convex hull
  of the rows of \code{M}.  So the question of interest becomes a 
  constrained optimization problem.
  
  Solving this problem relies on the package \code{lpSolve} to solve
  a linear program.  We may put the program in "standard form" by 
  writing \eqn{z=a-b}, where \eqn{a} and \eqn{b} are nonnegative
  vectors.  If we write \eqn{x=(a' b')'}, we obtain the linear program
  given by:
  
  Minimize \eqn{(-q' q')x} subject to \eqn{x'(L -L) \le 0} and
  \eqn{x \ge 0}.  One additional constraint arises because whenever any
  strictly negative value of \eqn{(-q' q')x} may be achieved, doubling
  \eqn{x} arbitrarily many times makes this value arbitrarily large
  in the negative direction, so no minimizer exists.  Therefore, we 
  add the constraint \eqn{(q' -q')x \le 1}.
  
  This function is used in the "stepping" algorithm of Hummel et al
  (2012).
}
\value{
  Logical, telling whether \code{p} is in the closed convex hull of the points
  in \code{M}.
}
\references{ 
\itemize{
    \item \url{http://www.cs.mcgill.ca/~fukuda/soft/polyfaq/node22.html}

    \item Hummel, R. M., Hunter, D. R., and Handcock, M. S. (2012),
    Improving Simulation-Based Algorithms for Fitting ERGMs,
    Journal of Computational and Graphical Statistics, 21: 920-939.
}
}

%\seealso{}
%\examples{}
%\keyword{  }% at least one, from doc/KEYWORDS
