%  File ergm/man/g4.Rd
%  Part of the statnet package, http://statnet.org
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
%  Copyright 2012 the statnet development team
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{g4}
\alias{g4}
\title{Goodreau's four node network as a ``network" object}
\description{
  This is an example thought of by Steve Goodreau. It is a directed
  network of four nodes and five ties stored as a \code{\link[network]{network}} object.

  It is interesting because the maximum likelihood estimator of the 
  model with out degree 3 in it exists, but the maximum psuedolikelihood
  estimator does not.
}
\usage{
 data(g4)
}
\source{
Steve Goodreau
}
\seealso{florentine, network, plot.network, ergm}
\examples{
data(g4)
summary(ergm(g4 ~ odegree(3), estimate="MPLE"))
summary(ergm(g4 ~ odegree(3), control=control.ergm(init=0)))
}
\keyword{datasets}

