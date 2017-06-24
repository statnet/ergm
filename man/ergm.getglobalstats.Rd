%  File man/ergm.getglobalstats.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2017 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm.getglobalstats}
\alias{ergm.getglobalstats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
internal function to return global statistics for a given network
}
\description{
The \code{ergm.getglobalstats} function is a low-level function not normally called by the user. It calculates and returns the global statistics for a given network and model.
}
\usage{
ergm.getglobalstats(nw, m, response = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nw}{
a \code{\link{network}} object
}
  \item{m}{
the model in use with network \code{nw}, as returned by \code{\link{ergm.getmodel}}
}
  \item{response}{
character name of an edge attribute to be used (for weighted ergm models)
}
}
\details{
 Calculates and returns the global statistics for a given network via \code{\link{ergm.Cprepare}} and \code{network_stats_wrapper.C} or \code{wt_network_stats_wrapper.C} if the model is weighted. It is called by \code{\link{summary.statistics.network}} which is generally the better way to access the functionality. 
}
\value{
returns a vector of the global statistics
}

\seealso{
\code{\link{summary.statistics.network}}
}
