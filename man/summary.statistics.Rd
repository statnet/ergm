\name{summary.statistics}
\alias{summary.statistics}
\alias{summary.statistics.matrix}
\alias{summary.statistics.network}
\alias{summary.statistics.default}
\alias{summary.statistics.formula}
\alias{summary.statistics.ergm}
\alias{summary.formula}
\alias{summary}
\title{Calculation of network or graph statistics}
\usage{
%summary(object, \dots)
%\method{summary.statistics}{default}(object,\dots, m=NULL)
\method{summary.statistics}{formula}(object, \dots, drop=FALSE, basis=NULL)
\method{summary.statistics}{ergm}(object, \dots, drop=FALSE, basis=NULL)
}
\description{
  Used to calculate the specified 
  statistics for an observed network if its argument is a formula for an 
  \code{\link{ergm}}.
  See \code{\link{ergm-terms}} for more information on the statistics that may
  be specified.
}
\arguments{
%  \item{formula}{formula; an \R \code{\link{formula}} object, of the form
%    \code{y ~ <model terms>},
%    where \code{y} is a \code{\link[network]{network}} object or a matrix that can be
%    coerced to a \code{\link[network]{network}}  object.  For the details on the possible
%    \code{<model terms>}, see \code{\link{ergm-terms}}.  To create a
%    \code{\link[network]{network}} object in \R, use the \code{network()} function,
%    then add nodal attributes to it using the \code{\%v\%}
%    operator if necessary.}
  \item{object}{an \R object. It is either an \R \code{\link{formula}} object
    (see above) or an 
    \code{\link{ergm}} model object.  
    In the latter case, \code{summary.statistics} is called for the
    \code{object$formula} object.
    In the former case, \code{object} is of the form
    \code{y ~ <model terms>},
    where \code{y} is a \code{\link[network]{network}} object or a matrix that can be
    coerced to a \code{\link[network]{network}}  object.  For the details on the possible
    \code{<model terms>}, see \code{\link{ergm-terms}}.  To create a
    \code{\link[network]{network}} object in \R, use the \code{network()} function,
    then add nodal attributes to it using the \code{\%v\%}
    operator if necessary.
    }
% \item{m}{A model component produced by \code{ergm.getmodel}. Only the 
%   \code{m$coef.names} component is used to set the names of the statistics.}
  \item{drop}{logical:  Should terms whose observed statistics are extreme
  among the set of all possible network statistics
  (which result in nonexistent MLEs) be dropped?}
  \item{basis}{An optional \code{\link[network]{network}} object relative to which
  the global statistics should be calculated.}
  \item{\dots}{further arguments passed to or used by methods.}
}

\details{
  If \code{object} is of class \code{\link{formula}}, then
  \code{\link[base]{summary}} may be used in lieu of \code{summary.statistics}
  because \code{summary.formula} calls the \code{summary.statistics} function.
 %If neither of those are given as the
 %first argument then a \code{ergm.Cprepare} object is expected. This last
 %option is meant for internal use.
 The function actually cumulates the change statistics when removing edges
 from the observed network one by one until the empty network results.  
 Since each model term has a prespecified value (zero by default) for the
 corresponding statistic(s) on an empty network, these change statistics
 give the absolute statistics on the original network.

 %More information can be found by looking at the documentation of
 %\code{\link{ergm}}.
}
\value{A vector of statistics measured on the network.
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
m <- as.matrix(flomarriage)
summary(m ~ edges)  # twice as large as it should be
summary(m ~ edges, directed=FALSE) # Now it's correct
}
\keyword{models}
