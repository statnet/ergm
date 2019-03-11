%  File man/as.edgelist.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, https://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  https://statnet.org/attribution
%
%  Copyright 2003-2019 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{as.edgelist}
\alias{as.edgelist.compressed}
\title{
Convert a network object into a numeric edgelist matrix
}
\description{
Constructs an edgelist in the format expected by ergm's internal functions

\strong{NOTE:} the \code{as.edgelist} functions have been moved to the
\code{network} package, and this help file may be removed in the future.  See \code{\link[network]{as.edgelist}}
}


\details{
Constructs an edgelist matrix from a network, sorted tails-major order, with tails first, and, for undirected networks, tail < head.

The \code{\link[network]{as.matrix.network}(nw, matrix.type="edgelist")} provides similar functionality but it does not enforce ordering..  


}



\note{
The \code{as.edgelist} functions have been moved to the \code{network} package.  See \code{\link[network]{as.edgelist}}
}



\seealso{
See also\code{\link[network]{as.edgelist}}, \code{\link[network]{as.matrix.network.edgelist}}
}
\examples{
\donttest{
data(faux.mesa.high)
as.edgelist(faux.mesa.high)
}
}
