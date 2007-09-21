\name{summary.gofobject}
\alias{summary.gofobject}
\alias{print.gofobject}
\title{Summaries the Goodness-of-Fit Diagnostics on a Exponential Family Random Graph Model}
\description{
  \code{\link{summary.gofobject}} summaries the diagnostics such as the 
  degree distribution, geodesic distances, shared partner distributions,
  and reachability for the goodness-of-fit of exponential family random graph
  models.  See \code{\link{ergm}} for more information on these models.
}
\usage{
\method{summary}{gofobject}(object, \dots)
}
\arguments{
  \item{object}{an object of class \code{gofobject}, 
       typically produced by the
       \code{\link{gof.ergm}} or \code{\link{gof.formula}} functions.
       See the documentation for these.}
  \item{\dots}{Additional arguments, to be passed to the plot function.}
}

\details{
  \code{\link{gof.ergm}} produces a sample of networks randomly drawn from the specified model.
  This function produces a print out the summary measures.
}
\value{none}

\seealso{gof.ergm, gof.formula, ergm, network, simulate.ergm}
\examples{
#
data(florentine)
#
# test the gof.ergm function
#
gest <- ergm(flomarriage ~ edges + kstar(2))
gest
summary(gest)

#
# Plot the probabilities first
#
gofflo <- gof(gest)
gofflo
summary(gofflo)
}
\keyword{graphs}
