\name{control.gof}
\alias{control.gof}
\alias{control.gof.formula}
\alias{control.gof.ergm}
\title{ Auxiliary for Controlling ERGM Goodness-of-Fit Evaluation  }
\description{
  Auxiliary function as user interface for fine-tuning ERGM
  Goodness-of-Fit Evaluation.
}
\usage{
control.gof.formula(prop.weights = "default", prop.args = NULL, drop =
TRUE, summarizestats = FALSE, maxchanges = 1e+06)

control.gof.ergm(prop.weights = NULL, prop.args = NULL, drop = TRUE, summarizestats = FALSE, maxchanges = 1e+06)
}
\details{
    This function is only used within a call to the \code{\link{gof}} function.
    See the \code{usage} section in \code{\link{gof}} for details.
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
    Default is \code{TRUE}.}
  \item{summarizestats}{logical; Print out a summary of 
    the sufficient statistics of the generated network. This is useful 
    as a diagnostic. Default is \code{FALSE}.}
  \item{maxchanges}{
    %FIXME: I am not sure this is accurate.
    %Maximum number of changes in dynamic network
    %simulation for which to allocate space. 
    Currently unused.}
}
\value{
  A list with arguments as components.
}
\seealso{ \code{\link{gof}}. The \code{\link{control.simulate}} 
function performs a 
similar function for
\code{\link{simulate.ergm}}; 
\code{\link{control.ergm}} performs a
similar function for \code{\link{ergm}}.   }
%\keyword{goodness of fit}
\keyword{ models }

