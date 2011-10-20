\name{control.simulate.stergm}
\alias{control.simulate.stergm}
\title{ Auxiliary for Controlling Separable Temporal ERGM Simulation }
\description{
  Auxiliary function as user interface for fine-tuning STERGM simulation.
}
\usage{
control.simulate.stergm(prop.weights = "default", prop.args = NULL, 
summarizestats = FALSE, 
output.network="network",
maxchanges = 1e+06, 
packagenames="ergm", 
parallel=0)
}
\details{
    This function is only used within a call to the \code{\link{simulate}} function.
    See the \code{usage} section in \code{\link{simulate.stergm}} for details.
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
  \item{summarizestats}{logical; Print out a summary of 
    the sufficient statistics of the generated network. This is useful 
    as a diagnostic. Default is \code{FALSE}.}
  \item{output.network}{character; class with which to store the output networks. The options are "network" (default) and "as.edgelist.compressed" (which saves space but only supports networks with vertex attributes}
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
\seealso{ \code{\link{simulate.stergm}}, \code{\link{simulate.formula}}. 
\code{\link{control.stergm}} performs a 
similar function for
\code{\link{stergm}}.   }
%\keyword{simulation}
\keyword{models}
