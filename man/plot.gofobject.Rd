%  File ergm/man/plot.gofobject.Rd
%  Part of the statnet package, http://statnet.org
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
%  Copyright 2012 the statnet development team
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{plot.gofobject}
\alias{plot.gofobject}
\title{Plot Goodness-of-Fit Diagnostics on a Exponential Family Random Graph Model}
\description{
  \code{\link{plot.gofobject}} plots diagnostics such as the 
  degree distribution, geodesic distances, shared partner distributions,
  and reachability for the goodness-of-fit of exponential family random graph
  models.  See \code{\link{ergm}} for more information on these models.
}
\usage{
\method{plot}{gofobject}(x, \dots,
         cex.axis=0.7, plotlogodds=FALSE,
         main = "Goodness-of-fit diagnostics",
         normalize.reachability=FALSE,
         verbose=FALSE)
}
\arguments{
  \item{x}{an object of class \code{gofobject}, 
       typically produced by the
       \code{\link{gof.ergm}} or \code{\link{gof.formula}} functions.
       See the documentation for these.}
  \item{cex.axis}{Character expansion of the axis labels relative to that for
   the plot.}
  \item{plotlogodds}{Plot the odds of a dyad having given
  characteristics (e.g., reachability, minimum geodesic
  distance, shared partners). This is an alternative to the
  probability of a dyad having the same property.}
  \item{main}{Title for the goodness-of-fit plots.}
  \item{normalize.reachability}{Should the reachability proportion be 
   normalized to make it more comparable with the other geodesic distance 
   proportions.}
  \item{verbose}{Provide verbose information on the progress of the
    plotting.}
  \item{\dots}{Additional arguments, to be passed to the plot function.}
}

\details{
  \code{\link{gof.ergm}} produces a sample of networks randomly drawn from the specified model.
  This function produces a plot of the summary measures.
}
\value{none}

\seealso{gof.ergm, gof.formula, ergm, network, simulate.ergm}
\examples{
\dontrun{
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
plot(gofflo)
#
# And now the odds 
#
plot(gofflo, plotlogodds=TRUE)
#
# Use the formula version
#
gof(flomarriage ~ edges + kstar(2), coef=c(-1.6339, 0.0049))
}
}
\keyword{graphs}
