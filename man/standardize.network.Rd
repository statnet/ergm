\name{standardize.network}
\alias{standardize.network}
\title{
Copy a network object enforcing ergm-appropriate guarantees about its internal representation
}
\description{
Create a copy of a \code{\link{network}} of interest with certain guarantees about its internal representation:
\itemize{
\item for every edge, the id of the 'tails' vertex is < id of the 'heads' vertex
\item no (tail,head) id pair has more than one edge ID associated with it (no multiplex edges)
}

}
\usage{
standardize.network(nw, preserve.eattr = TRUE)
}
\arguments{
  \item{nw}{
a \code{\link{network}} object to be copied
}
  \item{preserve.eattr}{
logical; should the edge attributes be preserved during the copying process (presumably slower)
}
}
\details{
This function is needed because the \code{\link{network}} object can support added non-directed edges in arbitrary order, as well as multiplex edges and hypergraphs (even if the network is not marked as such), which are not supported in the underlying ergm processes. Uses \code{\link{as.edgelist}} internally to make the conversion.
}
\value{
returns an ergm-appropriate network object. 
}

\note{
This function may be time expensive.  Also, the transformation is performed by deleting edges in initial network and re-adding them.
}

\examples{
test<-network.initialize(5,directed=FALSE)
test[2,1]<-1  #ergm wont like this
test$mel[[1]] # peek at internal representation

test2<-standardize.network(test) # enforce!
test2$mel[[2]]  # 1 and 2 have traded places

}

\keyword{internal}
