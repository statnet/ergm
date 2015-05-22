\name{as.edgelist}
\alias{as.edgelist}
\alias{as.edgelist.compressed}
\alias{as.edgelist.matrix}
\alias{as.edgelist.network}
\alias{as.edgelist.uncompressed}
\title{
Convert a network object into a numeric edgelist matrix
}
\description{
Constructs an edgelist in the format expected by ergm's internal functions
}
\usage{
\method{as.edgelist}{network}(x, attrname = NULL, as.sna.edgelist = FALSE, 
                            inverted = NULL, ...)
\method{as.edgelist}{matrix}(x, n, directed = TRUE, bipartite = FALSE, 
                            loops = FALSE, inverted = FALSE, ...)
}

\arguments{
  \item{x}{
a \code{network} object with additional class added indicating how it should be dispatched.
}
 \item{attrname}{optionally, the name of an edge attribute to use for edge values}
 \item{as.sna.edgelist}{logical; should the edgelist be returned in edgelist form expected by the sna package?}
 \item{inverted}{logical; value is passed to the 'inverted' flag on the edgelist returned}
 \item{n}{integer number of vertices in network, value passed to the 'n' flag on edgelist returned}
 \item{directed}{logical; is network directed, value passed to the 'directed' flag on edgelist returned}
 \item{bipartite}{logical or integer; is network bipartite, value passed to the 'bipartite' flag on edgelist returned}
 \item{loops}{logical; are self-loops allowed in network?, value passed to the 'loops' flag on edgelist returned}
  \item{\dots}{
additional arguments to other methods
}
}
\details{
Constructs an edgelist matrix from a network, sorted tails-major order, with tails first, and, for undirected networks, tail < head.

The \code{\link[network]{as.matrix.network}(nw, matrix.type="edgelist")} provides similar functionality but it does not enforce ordering..  


}
\value{
a two-column matrix of integers giving the tail (source) and head (target) indices of each edge.  The edgelist has additional attributes attached to it
\itemize{
\item \code{attr(,"n")} the number of vertices in the original network

\item \code{attr(,"directed")} logical, was the original network directed

\item \code{attr(,"bipartite")} was the original network bipartite

\item \code{attr(,"loops")} does the original network contain self-loops

\item \code{attr(,"inverted")}  indicates if the network has been inverted?
}
}


\note{
NOTE THAT THIS FUNCTION WILL PROBABLY BE TRANSFERED TO THE network PACKAGE IN THE NEAR FUTURE
}



\seealso{
See also \code{\link[network]{as.matrix.network.edgelist}}
}
\examples{
data(faux.mesa.high)
as.edgelist(faux.mesa.high)
}
