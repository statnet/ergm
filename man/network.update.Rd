\name{network.update}
\alias{network.update}
\title{Replaces the sociomatrix in a network object}
\usage{
network.update(nw, newmatrix)
}
\description{
  Replaces the sociomatrix in a network object with the sociomatrix specified by
\code{newmatrix}.
  See \code{\link{ergm}} for more information.
}
\arguments{
  \item{nw}{a \code{\link[network]{network}} object. See documentation for the \code{\link[network]{network}}
   package.}
  \item{newmatrix}{a matrix of zeros and one indicating the presence of a tie from i to j.}
}

\value{\code{\link{network.update}} returns a \code{\link[network]{network}} object.
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
network.update(flomarriage, rand.mat)
# Try this with an edgelist
rand.mat <- as.matrix.network.edgelist(flomarriage)[1:5,]
network.update(flomarriage, rand.mat)
}
\keyword{models}
