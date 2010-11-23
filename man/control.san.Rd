\name{control.san}
\alias{control.san}
\title{ Auxiliary for Controlling SAN }
\description{
  Auxiliary function as user interface for fine-tuning simulated annealing
  algorithm.
}
\usage{
control.san(prop.weights = "default", prop.args = NULL,
            drop = FALSE, maxedges=20000, maxchanges=1000000,
            packagenames="ergm", network.output="network",
            parallel = 0)
}
\details{
    This function is only used within a call to the \code{\link{san}} function.
    See the \code{usage} section in \code{\link{san}} for details.
}
\arguments{
  \item{prop.weights}{Specifies the method to allocate probabilities of
    being proposed to dyads. Defaults to \code{"default"}, which picks a
    reasonable default for the specified constraint. 
    Other possible values are
    \code{"TNT"}, \code{"random"}, and \code{"nonobserved"}, though not
    all values may be used
    with all possible constraints.}
  \item{prop.args}{An alternative, direct way of specifying additional 
    arguments to proposal.}
  \item{drop}{logical; Should the degenerate terms in the model be
    dropped from the fit?
    If statistics occur on the extreme of their range they
    correspond to infinite parameter estimates.
    Default is \code{FALSE}.}
  \item{maxedges}{maximum number of edges expected in network}
  \item{packagenames}{Names of packages in which changestatistics are found. Currently unused as \code{ergm} is presumed.}
  \item{network.output}{character: R class with which to output 
  networks. The options are "network" (default) and "edgelist.compressed" 
  (which saves space but only supports networks without vertex attributes)}
  \item{maxchanges}{
    %FIXME: I am not sure this is accurate.
    %Maximum number of changes in dynamic network
    %simulation for which to allocate space. 
    Currently unused}
  \item{parallel}{Number of threads in which to run the sampling. Currently unimplemented}
}
\value{
  A list with arguments as components.
}
\seealso{ \code{\link{san}}
}
%\keyword{simulation}
\keyword{models}
