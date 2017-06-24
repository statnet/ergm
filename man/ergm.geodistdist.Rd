%  File man/ergm.geodistdist.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2017 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm.geodistdist}
\alias{ergm.geodistdist}
\alias{ergm.geodesicmatrix}
\alias{ergm.geodesicmatrix.edgelist}
\alias{ergm.geodistn}
\alias{ergm.nodegeodesics}
\alias{ergm.pairgeodesic}
\title{
calculate geodesic distance distribution for a network or edgelist
}
\description{
\code{ergm.geodistdist} calculates geodesic  distance distribution for a given \code{\link{network}} and returns it as a vector.

\code{ergm.geodistn} calculates geodesic deistance distribution based on an input edgelist, and has very little error checking so should not normally be called by users. The C code requires the edgelist to be directed and sorted correctly.

}
\usage{
ergm.geodistdist(nw, directed = is.directed(nw))

ergm.geodistn(edgelist, n = max(edgelist), directed = FALSE)
}

\arguments{
  \item{nw}{
 \code{\link{network}} object over which distances should be calculated
}
  \item{directed}{
logical, should the network be treated as directed
}
  \item{edgelist}{an edgelist representation of a network as an mx2 matrix}
  \item{n}{integer, size of the network}
}
\details{
\code{ergm.geodistdist} is a network wrapper for \code{ergm.geodistn}, which calculates and returns the geodesic  distance distribution for a given network via full_geodesic_distribution.C
}
\value{
a vector \code{ans} with length equal to the size of the network where
\itemize{
     \item ans[i], i=1, ..., n-1 is the number of pairs of geodesic length i
     \item ans[n] is the number of pairs of geodesic length infinity.
     }
}



\seealso{
See also the sna package \code{\link[sna]{geodist}} function
}
\examples{
data(faux.mesa.high)
ergm.geodistdist(faux.mesa.high)
}
