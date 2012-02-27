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
control.gof.formula(nsim=100,
                              MCMC.burnin=1000,
                              MCMC.interval=1000,
                              MCMC.prop.weights="default",
                              MCMC.prop.args=list(),

                              MCMC.init.maxedges=20000,
                              MCMC.packagenames="ergm",

                              MCMC.runtime.traceplot=FALSE,          
                              summarizestats=FALSE,
                              network.output="network",

                              seed=NULL,
                              parallel=0,
                              parallel.type=NULL,
                              parallel.version.check=TRUE)

control.gof.ergm(nsim=100,
                 MCMC.burnin=NULL,
                 MCMC.interval=NULL,
                 MCMC.prop.weights=NULL,
                 MCMC.prop.args=NULL,
                 MCMC.init.maxedges=NULL,
                 MCMC.packagenames=NULL,

                 MCMC.runtime.traceplot=FALSE,
                 summarizestats=FALSE,
                 network.output="network",

                 seed=NULL,
                 parallel=0,
                 parallel.type=NULL,
                 parallel.version.check=TRUE)
}
\details{
    This function is only used within a call to the \code{\link{gof}} function.
    See the \code{usage} section in \code{\link{gof}} for details.
}
\arguments{
  \item{MCMC.prop.weights}{Specifies the method to allocate probabilities of
    being proposed to dyads. For the \code{\link{simulate.formula}}
    variant, defaults to \code{"default"}, which picks a
    reasonable default for the specified constraint. For
    \code{\link{simulate.ergm}} variant, defaults to \code{NULL}, to
    reuse the weights with which the given \code{\link{ergm.object}} was
    fitted. Other possible values are
    \code{"TNT"}, \code{"random"}, and \code{"nonobserved"}, though not
    all values may be used
    with all possible constraints.}
  \item{MCMC.prop.args}{An alternative, direct way of specifying additional arguments to proposal.}
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

