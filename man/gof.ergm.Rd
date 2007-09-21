\name{gof.ergm}
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
         proposaltype="randomtoggle",
         seed=NULL, drop=TRUE, boundDeg = NULL, summarizestats=FALSE,
         verbose=FALSE)
\method{gof}{ergm}(object, \dots,
         nsim=100,
         GOF=~degree+espartners+distance,
         burnin=10000, interval=1000,
         proposaltype="randomtoggle",
         seed=NULL, drop=TRUE, summarizestats=FALSE,
         theta0=NULL, verbose=FALSE)
}
\arguments{
 \item{object}{an R object. Either a formula or an ergm object. 
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
 \item{proposaltype}{The type of Metropolis-Hastings proposal used.
   The default is \code{"randomtoggle"} to propose a dyad at random
   to toggle (i.e., change the value of the dyad to its complement).
   The options are 
    \describe{
     \item{\code{TNT}}{("tie/no-tie") proposes a dyad with
        at least 50\% chance of toggling an existing   
        edge, as opposed to a simple random toggle that rarely does so
        in sparse networks.}
     \item{\code{ConstantEdges}}{proposes a pair of dyads with one having a tie
        and the other not having one.}
     \item{\code{CondDegree}}{proposes a pair of dyads in such a way
        to retain the individual degrees of the nodes.}
     \item{\code{CondDegreeDist}}{proposes either a hexad or tetrad of dyads
        to toggle in a way to retain the degree distribution of the 
        network (but not necessarily the individual degrees).}
     \item{\code{CondOutDegreeDist}}{proposes either a hexad or tetrad of dyads
        to toggle in a way to retain the out degree distribution of the 
        network (but not necessarily the individual degrees).}
     \item{\code{CondInDegreeDist}}{proposes either a hexad or tetrad of dyads
        to toggle in a way to retain the in degree distribution of the 
        network (but not necessarily the individual degrees).}
     \item{\code{RandomNode}}{to propose all dyads associated with a randomly
        sampled node.}
     \item{\code{TwoRandomToggles}}{proposes a pair of randomly chosen
        dyads.}
    }
  The option \code{"CondDegree"} is forced when
  \code{condAllDegExact=TRUE}.
  This is essential for some models or \code{boundDeg} constraints.}
% \item{multiplicity}{The number of replications of proposal types
%  grouped together to form a single proposal.}
  \item{drop}{Should the degenerate terms in the model be   
   dropped from the fit?
   If statistics occur on the extreme of their range they
   correspond to infinite parameter estimates.}
  \item{boundDeg}{A list comprising of up to five matrices and one or more
    arguments, which are used to constrain the number of in or out
    degrees that a node can have. By default no constraining is applied.
    This is explained more fully in the \code{\link{ergm}} function.}
  \item{summarizestats}{Print out a summary of the sufficient statistics of
   the generated graph. This is useful as a diagnostic.}
  \item{seed}{Random number integer seed.
   The default is \code{sample(10000000, size=1)}.}
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
#
data(florentine)
#
# test the gof.ergm function
#
gest <- ergm(flomarriage ~ edges + kstar(2))
gest
summary(gest)

#
# Plot the probabilities first
#
gofflo <- gof(gest)
gofflo
#
# Place all three on the same page
# with nice margins
#
par(mfrow=c(1,3))
par(oma=c(0.5,2,1,0.5))
#
plot(gofflo)
#
# And now the odds 
#
plot(gofflo, plotodds=TRUE)
#
# Use the formula version
#
plot(gof(flomarriage ~ edges + kstar(2), theta0=c(-1.6339, 0.0049)))
}
\keyword{models}
