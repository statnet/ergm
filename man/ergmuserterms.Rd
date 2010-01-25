\name{ergmuserterms-package}
\alias{ergmuserterms}
\alias{ergmuserterms}
\docType{package}
\title{Add Statistics Terms for the 'ergm' Package}
\description{
  The \code{\link{ergm}} package is capable of fitting a wide range of 
  exponential random network models, in which
  the probability of a given network, \eqn{y}, on a set of nodes is 
  \eqn{\exp(\theta{\cdot}g(y))/c(\theta)}, where 
  \eqn{g(y)} is a vector of network statistics,
  \eqn{\theta} is a parameter vector of the same length and \eqn{c(\theta)} is the
  normalizing constant for the distribution.
  The \code{\link{ergm}} function fits these models when they are expressed
  via an \R \code{\link{formula}} object, of the form
 \code{y ~ <model terms>},
  where \code{y} is a network object or a matrix that can be
  coerced to a network object.  
  To create a
  network object in \R, use the \code{network()} function,
  then add nodal attributes to it using the \code{\%v\%}
  operator if necessary.

  The \code{\link{ergm}} package contains a wide range of terms.
  For the details on the possible
  \code{<model terms>}, see \code{\link{ergm-terms}}.  

  This package can be modified by users to add user-defined terms to \code{\link{ergm}}
  models.
  The terms can be used throughout the \code{\link{ergm}} package
  and behave identically to the supplied terms.
}
%\usage{
%terms.ergm()
%}
%\arguments{
%\item{\dots}{Specific to the model term.}
%}

\details{
    The \code{ergmuserterms} package is available from the \code{statnet}
    website (\url{http://statnetproject.org}).

    The code contains some simple examples and templates.
    These include:
    \itemize{
      \item{\code{m2star}}{\emph{Mixed 2-stars, a.k.a. 2-paths.} 
      This option can only be 
      specified with a directed network; for undirected graphs 
      see \code{kstar(2)}. This 
      option adds one statistic to the model, equal to the 
      number of mixed-2-stars in the network, defined as a 
      pair of edges \eqn{\{(i{\rightarrow}j), (j{\rightarrow}k)\}.}
      }

      \item{\code{testme}}{\emph{A clone of Edges.} 
      This is included for purposes of an example.
      This option adds one graph statistic
      equal to the number of edges in the graph.  For undirected
      graphs, \code{edges} is isomorphic to \code{kstar(1)}; for
      directed networks, \code{edges} is isomorphic to both \code{ostar(1)}
      and \code{istar(1)}.  }
    }
  In the implementation of \code{\link{ergm}}, the model is initialized in \R,
  then all the model information is passed to a \code{C} program that generates
  the sample of graph statistics using MCMC.  This sample is then returned
  to \R, which then approximates the MLE.
}
%\value{\code{\link{ergm}} returns an object of class 'ergm' that is a list.
%}
\seealso{ergm, network, ergm-terms}
\examples{
\dontrun{
library(ergmuserterms)
data(sampson)
monk.fit <- ergm(samplike~m2star)
summary(monk.fit)

monk.fit <- ergm(samplike ~ m2star + testme)
summary(monk.fit)
}
}
\keyword{graphs}
