\name{ergm.simulate.control}
\alias{ergm.simulate.control}
\title{ Auxiliary for Controlling MCMC Simulation}
\description{
  Auxiliary function as user interface for fine-tuning MCMC simulation of an ergm.
}
\usage{
ergm.simulate.control(prop.weights = "default", prop.args = NULL, boundDeg = NULL, 
                      drop = FALSE, summarizestats = FALSE, maxchanges = 1e+06)
}
\arguments{
  \item{prop.weights}{ Specifies the method to allocate probabilities of
    being proposed to dyads. Defaults to \code{"default"}, which picks a
    reasonable default for the specified constraint. Possible values are
    "TNT", "random", and "nonobserved", though not all values may be used
    with all possible constraints (in the \code{\link{ergm}} function)
    }
  \item{prop.args}{ An alternative, direct way of specifying additional arguments to proposal. }
  \item{\code{boundDeg}}{list; comprising up to five matrices and 
    one or more arguments, which are used to constrain the number of in 
    or out degrees that a node can have. By default no constraining 
    is applied.  This is explained more fully under the 
    \code{Conditioning within models} section of documentation for \code{\link{ergm}}. Defaults to no conditioning.}
  \item{\code{drop}}{logical; Should the degenerate terms in the model be
    dropped from the fit?
    If statistics occur on the extreme of their range they
    correspond to infinite parameter estimates.
    Default is \code{TRUE}.}
  \item{summarizestats}{ DESCRIBE SUMMARIZESTATS HERE}
  \item{maxchanges}{ DESCRIBE MAXCHANGES HERE  }
}
\value{
  A list with arguments as components.
}
\seealso{ \code{\link{ergm}}, \code{\link{simulate.ergm}}, \code{\link{ergm.control}},
  \code{\link{glm.control}} performs a
  similar function for \code{\link{glm}} }

\keyword{ models }

