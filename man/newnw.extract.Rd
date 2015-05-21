\name{newnw.extract}
\alias{newnw.extract}
\title{
Internal function to create a new network from the ergm MCMC sample output
}
\description{
An internal function to generate a new \code{\link{network}} object using the output (lists of toggled heads and tail vertices) from an ergm MCMC or SAN process. 
}
\usage{
newnw.extract(oldnw, z, output = "network", response = NULL)
}
\arguments{
  \item{oldnw}{
a network object (presumably input to the ergm process) from which the network- and vertex-level attributes will be copied
}
  \item{z}{
a list having either a component named \code{newedgelist} or two components \code{newtails} and \code{newheads} containing the ids of the head and tails vertices of the edges.  Optionall \code{newweights} containing edgewights.
}
  \item{output}{
passed to \code{\link{network.update}}, which claims not to use it
}
  \item{response}{
optional character string giving the name of the edge attribute where the edge values (weight/count) should be stored.
}
}

\value{
a \code{\link{network}} object with properties copied from \code{oldnw} and edges corresponding to the lists of tails and head vertex ids in \code{z}
}

\note{
This is an internal ergm function, it most cases with edgelists to be converted to networks it will probably be simpler to use \code{\link{network.edgelist}}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{network.edgelist}}, \code{\link{network.update}}
}
