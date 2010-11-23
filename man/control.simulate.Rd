\name{control.simulate}
\alias{control.simulate}
\alias{control.simulate.formula}
\alias{control.simulate.ergm}
\title{ Auxiliary for Controlling ERGM Simulation }
\description{
  Auxiliary function as user interface for fine-tuning ERGM simulation.
}
\usage{
control.simulate(prop.weights = "default", prop.args = NULL,
drop = FALSE, summarizestats = FALSE, maxedges=20000,
maxchanges = 1e+06, 
packagenames="ergm", 
network.output="network",
parallel=0)

control.simulate.formula(prop.weights = "default", prop.args = NULL,
drop = FALSE, summarizestats = FALSE, maxedges=20000,
maxchanges = 1e+06, 
packagenames="ergm", 
network.output="network",
parallel=0)

control.simulate.ergm(prop.weights = "default", prop.args = NULL, drop = FALSE, 
summarizestats = FALSE, 
maxchanges = 1e+06, maxedges = 20000,
packagenames="ergm", 
network.output="network",
parallel=0)
}
\details{
    This function is only used within a call to the \code{\link{simulate}} function.
    See the \code{usage} section in \code{\link{simulate.ergm}} for details.
}
\arguments{
  \item{prop.weights}{Specifies the method to allocate probabilities of
    being proposed to dyads. For the \code{\link{simulate.formula}}
    variant, defaults to \code{"default"}, which picks a
    reasonable default for the specified constraint. For
    \code{\link{simulate.ergm}} variant, defaults to \code{NULL}, to
    reuse the weights with which the given \code{\link{ergm.object}} was
    fitted. Other possible values are
    \code{"TNT"}, \code{"random"}, and \code{"nonobserved"}, though not
    all values may be used
    with all possible constraints.}
  \item{prop.args}{An alternative, direct way of specifying additional arguments to proposal.}
  \item{drop}{logical; Should the degenerate terms in the model be
    dropped from the fit?
    If statistics occur on the extreme of their range they
    correspond to infinite parameter estimates.
    Default is \code{FALSE}.}
  \item{summarizestats}{logical; Print out a summary of 
    the sufficient statistics of the generated network. This is useful 
    as a diagnostic. Default is \code{FALSE}.}
  \item{maxedges}{maximum number of edges expected in network}
  \item{network.output}{character; class with which to store the output 
  networks. The options are "network" (default) and "edgelist.compressed" 
  (which saves space but only supports networks without vertex attributes)}
  \item{maxchanges}{
    %FIXME: I am not sure this is accurate.
    %Maximum number of changes in dynamic network
    %simulation for which to allocate space. 
    Currently unused}
  \item{packagenames}{Names of packages in which changestatistics are found.}
  \item{parallel}{ Number of threads in which to run the sampling. }
}
\value{
  A list with arguments as components.
}
\seealso{ \code{\link{simulate.ergm}}, \code{\link{simulate.formula}}. 
\code{\link{control.ergm}} performs a 
similar function for
\code{\link{ergm}};
\code{\link{control.gof}} performs a
similar function for \code{\link{gof}}.   }
%\keyword{simulation}
\keyword{models}
