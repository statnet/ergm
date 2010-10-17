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
is.inCH(p, M)
}
\arguments{
  \item{p}{A \eqn{d}-dimensional vector}
  \item{M}{An \eqn{r} by \eqn{d} matrix.  Each row of \code{M} is a 
  \eqn{d}-dimensional vector.}
}
\details{
  This function depends on the package \code{Rglpk} to solve a constrained
  linear optimization problem in order to determine an answer.
  The question of whether \code{p} is in a closed convex set \code{S} 
  may be formulated as the question of whether there exists a separating 
  hyperplane between \code{p} and \code{S}, which may in turn be formulated
  as the question of whether the maximum possible value of a linear function,
  subject to constraints, has a strictly positive solution.
  
  Note that the answer given could be incorrect simply due to rounding error if the
  true maximum is close to zero.  For this reason, the 
  package \code{rcdd}, which produces exact rational-number
  solutions to linear programs, could be used instead of \code{Rglpk}.  
  However, this approach would require more computing and would therefore
  be slower.
}
\value{
  Logical, telling whether \code{p} is in the closed convex hull of the points
  in \code{M}.
}
\references{ 
\url{http://www.cs.mcgill.ca/~fukuda/soft/polyfaq/node22.html}
}

%\seealso{}
%\examples{}
%\keyword{  }% at least one, from doc/KEYWORDS
