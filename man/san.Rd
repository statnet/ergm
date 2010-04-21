\name{san}
\alias{san.ergm}
\alias{san.formula}
\alias{san.default}
\alias{san}
\title{Use Simulated Annealing to attempt to match a network to a vector
of mean statistics}
\usage{
\method{san}{formula}(object, nsim=1, seed=NULL, theta0=NULL,
                      tau=1, invcov=NULL,
                      burnin=10000, interval=10000,
                      meanstats=NULL,
                      basis=NULL,
                      sequential=TRUE,
                      constraints = ~.,
                      control = control.san(),
                      verbose=FALSE, \dots)
\method{san}{ergm}(object, nsim=1, seed=NULL, theta0=object$coef,
                      burnin=10000, interval=10000, 
                      meanstats=NULL,
                      basis=NULL,
                      sequential=TRUE,
                      constraints = NULL,
                      control = control.san(),
                      verbose=FALSE, \dots)
}
\description{
  This function attempts to find a network or networks
  whose statistics match those passed
  in via the \code{meanstats} vector.
}
\arguments{
 \item{object}{an \R object. Either
    a \code{\link{formula}} or an \code{\link{ergm}} object.
    The \code{\link{formula}} should be of the form \code{y ~ <model terms>},
    where \code{y} is a network object or a matrix that can be
    coerced to a \code{\link[network]{network}} object.  For the details 
    on the possible
    \code{<model terms>}, see \code{\link{ergm-terms}}.  To create a
    \code{\link[network]{network}} object in \R, use the 
    \code{network()} function,
    then add nodal attributes to it using the \code{\%v\%}
    operator if necessary.}
 \item{nsim}{Number of desired networks.}
 \item{seed}{Random number integer seed.}
 \item{theta0}{Parameter values used for MCMC simulations.}
 \item{tau}{Currently unused.}
 \item{invcov}{Initial inverse covariance matrix used to calculate
    Mahalanobis distance in determining how far a proposed MCMC move is from
    the \code{meanstats} vector.  If NULL, taken to be the covariance
    matrix returned when fitting the MPLE (if \code{theta0==NULL}) or the
    identity (otherwise).}
 \item{burnin}{Number of MCMC steps prior to recording first vector of 
    network statistics.}
 \item{interval}{Number of MCMC steps between recordings of network statistics}
 \item{meanstats}{A vector of the same length as the number of terms implied by
    the formula, which is either \code{object} itself in the case of
    \code{san.formula} or \code{object$formula} in the case of \code{san.ergm}.}
 \item{basis}{If not NULL, a network that forms the beginning of the
    Markov chain.  If NULL, this is taken to be the network named in the
    formula.}
 \item{sequential}{Logical:  Should the returned draws use the prior draw 
   as the starting network or always use the initially passed network?
   }
 \item{constraints}{A one-sided formula specifying one or more constraints
    on the support of the distribution of the networks being
    simulated. See the documentation for a similar argument for
    \code{\link{ergm}} for more information. For
    \code{simulate.formula}, defaults to no constraints. For
    \code{simulate.ergm}, defaults to using the same constraints as
    those with which \code{object} was fitted.}
 \item{control}{A list of control parameters for algorithm
   tuning. Constructed using \code{\link{control.san}}.}
 \item{verbose}{If this is \code{TRUE}, we will print out more information as 
  we run the program, including (currently) some goodness of fit 
  statistics.}
 \item{\dots}{Further arguments passed to or used by methods.}
}

%\details{}

\value{
  A network or list of networks that hopefully have network statistics close
  to the \code{meanstats} vector.
}
%\seealso{}
%\examples{}
\keyword{models}
