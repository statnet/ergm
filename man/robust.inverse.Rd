\name{robust.inverse}
\alias{robust.inverse}

\title{
Find the inverse of a matrix
}
\description{
The \code{robust.inverse} function attempts to return the inverse of a matrix H
 by either direct means (via R's \code{\link{solve}}) or via computation from the
 singular decomposition of H
}
\usage{
robust.inverse(H, tol = sqrt(.Machine$double.eps))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{H}{
a matrix, presumably a Hessian
}
  \item{tol}{
the tolerance, used for determining which of the singular values
       from the decompostion are postive; default=sqrt(.Machine$double.eps)
}
}

\value{
returns one of:
\itemize{
  \item{ H : the original matrix, if H could not be inverted directly and the
       singular decomposition of H incurred an error}
   \item{0 : a p by q matrix of 0's, where H is q by p, if none of the singular
       values from the decomposition are positive}
   \item{iH: the inverse of H, computed directly if possible, otherwise computed
       from the components of the singular decomposition}
       }
}

\seealso{
\code{\link{svd}}, \code{\link{solve}}
}
