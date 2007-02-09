\name{summary.statistics}
\alias{summary.statistics}
\alias{summary.statistics.network}
\alias{summary.statistics.formula}
\alias{summary.statistics.ergm}
\title{Calculation of network or graph statistics}
\usage{
%summary(object, \dots)
%\method{summary.statistics}{default}(object,\dots, m=NULL)
\method{summary.statistics}{network}(object, \dots, basis=NULL)
\method{summary.statistics}{ergm}(object, \dots, basis=NULL)
}
\description{
  \code{\link[base]{summary}} is used to calculate the specified 
  statistics for an observed network. The statistics are specified 
  by a formula.
  See \code{\link{ergm}} for more information on these models.
}
\arguments{
  \item{object}{an \R object. It is either an \R formula object, of the form
    \code{y ~ <model terms>},
     where \code{y} is a network object or a matrix that can be
     coerced to a network object.  For the details on the possible
     \code{<model terms>}, see \code{\link{ergm-terms}}.  To create a
     network object in \R, use the \code{network()} function,
     then add nodal attributes to it using the \code{\%v\%}
     operator if necessary.
     The statistics are calculated for the network \code{y}.
    Otherwise, it is a \code{\link{ergm}} model object, and the statistics are
    computed for the model \code{\link{formula}}.}
% \item{m}{A model component produced by \code{ergm.getmodel}. Only the 
%   \code{m$coef.names} component is used to set the names of the statistics.}
  \item{basis}{An optional \code{\link[network]{network}} object to calculate the global
    statistics relative to.}
  \item{\dots}{further arguments passed to or used by methods.}
}

\details{
  The statistics are calculated for the specified model.  The
  model is specified by the first argument of the function.  If
  the first argument is a formula then this defines the model.  If the
  first argument is the output of a call to \code{\link{ergm}} then the model
  used for that call is the one used.
 %If neither of those are given as the
 %first argument then a \code{ergm.Cprepare} object is expected. This last
 %option is meant for internal use.

  More information can be found by looking at the documentation of
  \code{\link{ergm}}.
}
\value{\code{\link{summary.statistics}} returns a vector of statistics corresponding to the model.
}
\seealso{ergm, network, ergm-terms}
\examples{
#
# Lets look at the Florentine marriage data
#
data(florentine)
#
# test the summary.statistics function
#
summary(flomarriage ~ edges + kstar(2))
}
\keyword{models}
