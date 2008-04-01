\name{gof}
\alias{gof.ergm}
\alias{gof.default}
\alias{gof.formula}
\alias{gof.ergm}
\alias{gof}
\title{Conduct Goodness-of-Fit Diagnostics on a Exponential Family Random Graph Model}
\description{
  \code{\link{gof}} calculates \eqn{p}-values for geodesic
  distance, degree, and reachability summaries to
  diagnose the goodness-of-fit of exponential family random graph
  models.  See \code{\link{ergm}} for more information on these models.
}
\usage{
\method{gof}{default}(object,\dots)
\method{gof}{formula}(formula, \dots, theta0=NULL,
         nsim=100, burnin=10000, interval=1000,
         GOF=~degree+espartners+distance,
         constraints=~.,
         control=control.gof.formula(),
         seed=NULL,
         verbose=FALSE)
\method{gof}{ergm}(object, \dots,
         nsim=100,
         GOF=~degree+espartners+distance,
         burnin=10000, interval=1000,
         constraints=NULL,
         control=control.gof.ergm(),
         seed=NULL, 
         theta0=NULL, verbose=FALSE)
}
\arguments{
 \item{object}{an \R object. Either a formula or an \code{\link{ergm}} object. 
   See documentation for \code{\link{ergm}}.}
 \item{formula}{formula; An \R formula object, of the form
 \code{y ~ <model terms>},
  where \code{y} is a network object or a matrix that can be
  coerced to a network object. This specifies the model to simulate from.
  For the details on the possible
  \code{<model terms>}, see \code{\link{ergm-terms}}.  To create a
  network object in \R, use the \code{network()} function,
  then add nodal attributes to it using the \code{\%v\%}
  operator if necessary.}
 \item{theta0}{
    When given either
    a formula or an object of class ergm, theta0 are the parameters
    from which the sample is drawn. By default set to a vector of 0.}
%\item{theta0missing}{
%   Optional parameters from which the missing data
%   sample is drawn. By default set to \code{theta0}.}
  \item{nsim}{The number of simulations to use for the MCMC \eqn{p}-values.
   This is the size of the sample of graphs to be randomly drawn
   from the distribution specified by the object on the set of all graphs.
   }
  \item{GOF}{formula; an \R formula object, of the form
   \code{~ <model terms>} specifying the
   statistics to use to diagnosis the goodness-of-fit of the model. 
   They do not need to be in the model formula specified in
   \code{formula}, and typically are not.
   Examples are the degree distribution ("degree"), 
   minimum geodesic distances ("dist"), and shared partner distributions
   ("espartners" and "dspartners").
   For the details on the possible
   \code{<model terms>}, see \code{\link{ergm-terms}}.}
  \item{burnin}{Number of proposed edge toggles before any
    MCMC sampling is done. If the model is correct this can be 0.
    Currently, there is no support for any check
    of the Markov chain mixing, so burnin should be set to a
    fairly large number.}
  \item{interval}{Number of proposed edge toggles between sampled statistics.
    The program prints a warning if too few proposed toggles are
    being accepted (if the number of proposed toggles between
    sampled observations ever equals an integral multiple of
    100*(1+the number of toggles accepted)).}
  \item{constraints}{A one-sided formula specifying one or more constraints
    on the support of the distribution of the networks being
    modeled. See the help for similarly-named argument in
    \code{\link{ergm}} for more information. For
    \code{gof.formula}, defaults to unconstrained. For \code{gof.ergm},
    defaults to the constraints with which \code{object} was fitted.}
  \item{control}{A list to control parameters, constructed using
    \code{\link{control.gof.formula}} or \code{\link{control.gof.ergm}}
    (which have different defaults).}
  \item{seed}{integer; random number integer seed. Defaults to \code{NULL} to
    use whatever the state of the random number generater is at the time
    of the call.}
  \item{verbose}{Provide verbose information on the progress of the
    simulation.}
  \item{\dots}{Additional arguments, to be passed to lower-level functions
  in the future.}
}

\details{
  A sample of graphs is randomly drawn from the specified model.
  The first argument is typically 
  the output of a call to \code{\link{ergm}} and the model
  used for that call is the one fit.

  A plot of the summary measures is plotted.
  More information can be found by looking at the documentation of
  \code{\link{ergm}}.
}
\value{\code{\link{gof}}, \code{\link{gof.ergm}}, and \code{\link{gof.formula}}
         return an object of class \code{gofobject}.
       This is a list of the tables of statistics and \eqn{p}-values.
       This is typically plotted using \code{\link{plot.gofobject}}.
      }

\seealso{ergm, network, simulate.ergm, summary.ergm, plot.gofobject}
\examples{
data(florentine)
gest <- ergm(flomarriage ~ edges + kstar(2))
gest
summary(gest)

# test the gof.ergm function
gofflo <- gof(gest)
gofflo

# Plot all three on the same page
# with nice margins
par(mfrow=c(1,3))
par(oma=c(0.5,2,1,0.5))
plot(gofflo)

# And now the log-odds
plot(gofflo, plotlogodds=TRUE)

# Use the formula version of gof
gofflo2 <-gof(flomarriage ~ edges + kstar(2), theta0=c(-1.6339, 0.0049))
plot(gofflo2)
}
\keyword{models}
