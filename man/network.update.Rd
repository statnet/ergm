%  File man/network.update.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2017 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{network.update}
\alias{network.update}
\title{Replaces the sociomatrix in a network object}
\usage{
network.update(nw, newmatrix, matrix.type=NULL, output="network",
              ignore.nattr = c("bipartite", "directed", "hyper", 
              "loops", "mnext", "multiple", "n"), ignore.vattr = c())
}
\description{
  Replaces the edges in a network object with the edges corresponding to the sociomatrix specified by
\code{newmatrix}.
  See \code{\link{ergm}} for more information.
}
\arguments{
  \item{nw}{a \code{\link[network]{network}} object. See documentation for the \code{\link[network]{network}}
   package.}
  \item{newmatrix}{Either an adjacency matrix (a matrix of zeros and ones indicating the presence of a 
  tie from i to j) or an edgelist (a two-column matrix listing origin and destination node numbers
  for each edge; note that in an undirected matrix, the first column should be the smaller
  of the two numbers).}
  \item{matrix.type}{One of "adjacency" or "edgelist" telling which type of
  matrix \code{newmatrix} is.  Default is to use the 
  \code{\link[network]{which.matrix.type}} function.}
  \item{output}{Currently unused.}
  \item{ignore.nattr}{character vector of the names of network-level attributes to ignore when updating network objects (defaults to standard network properties)}
  \item{ignore.vattr}{character vector of the names of vertex-level attributes to ignore when updating network objects}
}

\value{\code{\link{network.update}} returns a new \code{\link[network]{network}} object with the edges specified by \code{newmatrix} and network and vertex attributes copied from the input network \code{nw}. Input network is not modified. 
}
\seealso{ergm, network}
\examples{
#
data(florentine)
#
# test the network.update function
#
# Create a Bernoulli network
rand.net <- network(network.size(flomarriage))
# store the sociomatrix 
rand.mat <- rand.net[,]
# Update the network
network.update(flomarriage, rand.mat, matrix.type="adjacency")
# Try this with an edgelist
rand.mat <- as.matrix.network.edgelist(flomarriage)[1:5,]
network.update(flomarriage, rand.mat, matrix.type="edgelist")
}
\keyword{models}
