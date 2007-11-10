\name{simulate.ergm.control}
\alias{simulate.ergm.control}
\title{ Auxiliary for Controlling MCMC Simulation}
\description{
  Auxiliary function as user interface for fine-tuning MCMC simulation of an ergm.
}
\usage{
simulate.ergm.control(prop.weights = "default", prop.args = NULL,
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

