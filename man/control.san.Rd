\name{control.san}
\alias{control.san}
\title{ Auxiliary for Controlling SAN }
\description{
  Auxiliary function as user interface for fine-tuning simulated annealing
  algorithm.
}
\usage{
control.san(prop.weights = "default", prop.args = NULL,
            drop = FALSE, summarizestats = FALSE, 
            maxchanges = 1000000)
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
  \item{maxchanges}{
    %FIXME: I am not sure this is accurate.
    %Maximum number of changes in dynamic network
    %simulation for which to allocate space. 
    Currently unused}
}
\value{
  A list with arguments as components.
}
\seealso{ \code{\link{san}}
}
%\keyword{simulation}
\keyword{models}
