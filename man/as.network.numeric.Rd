% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.network.numeric.R
\name{as.network.numeric}
\alias{as.network.numeric}
\title{Create a Simple Random network of a Given Size}
\usage{
\method{as.network}{numeric}(
  x,
  directed = TRUE,
  hyper = FALSE,
  loops = FALSE,
  multiple = FALSE,
  bipartite = FALSE,
  ignore.eval = TRUE,
  names.eval = NULL,
  edge.check = FALSE,
  density = NULL,
  init = NULL,
  numedges = NULL,
  ...
)
}
\arguments{
\item{x}{count; the number of nodes in the network}

\item{directed}{logical; should edges be interpreted as directed?}

\item{hyper}{logical; are hyperedges allowed? Currently ignored.}

\item{loops}{logical; should loops be allowed? Currently ignored.}

\item{multiple}{logical; are multiplex edges allowed? Currently ignored.}

\item{bipartite}{count; should the network be interpreted as bipartite? If
present (i.e., non-NULL) it is the count of the number of actors in the
bipartite network. In this case, the number of nodes is equal to the number
of actors plus the number of events (with all actors preceding all events).
The edges are then interpreted as nondirected.}

\item{ignore.eval}{logical; ignore edge values? Currently ignored.}

\item{names.eval}{optionally, the name of the attribute in which edge values
should be stored. Currently ignored.}

\item{edge.check}{logical; perform consistency checks on new edges?}

\item{density}{numeric; the probability of a tie for Bernoulli networks. If
neither density nor \code{init} is given, it defaults to the number of nodes
divided by the number of dyads (so the expected number of ties is the same
as the number of nodes.)}

\item{init}{numeric; the log-odds of a tie for Bernoulli networks.  It is
only used if density is not specified.}

\item{numedges}{count; if present, sample the Bernoulli network conditional
on this number of edges (rather than independently with the specified
probability).}

\item{...}{additional arguments}
}
\value{
An object of class \code{\link[network:network]{network}}
}
\description{
\code{\link[=as.network.numeric]{as.network.numeric()}} creates a random Bernoulli network of the
given size as an object of class \code{\link[network:network]{network}}.
}
\details{
The network will not have vertex, edge or network attributes.  These
can be added with operators such as \code{\%v\%}, \code{\%n\%}, \code{\%e\%}.
}
\examples{
# Draw a random directed network with 25 nodes
g <- network(25)

# Draw a random undirected network with density 0.1
g <- network(25, directed=FALSE, density=0.1)

# Draw a random bipartite network with 4 actors and 6 events and density 0.1
g <- network(10, bipartite=4, directed=FALSE, density=0.1)

# Draw a random directed network with 25 nodes and 50 edges
g <- network(25, numedges=50)
}
\references{
Butts, C.T.  2002.  ``Memory Structures for Relational Data in
R: Classes and Interfaces'' Working Paper.
}
\seealso{
\code{\link[network:network]{network}}
}
\keyword{classes}
\keyword{graphs}
