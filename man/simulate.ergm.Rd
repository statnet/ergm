\name{simulate.ergm}
\alias{simulate.ergm}
\alias{simulate.formula}
\alias{simulate.formula.ergm}
\title{Draw from the distribution of an Exponential Family Random Graph Model}
\usage{
\method{simulate}{formula}(object, nsim=1, seed=NULL, coef,
                          burnin=1000, interval=1000,
                          basis=NULL,
                          statsonly=FALSE,
                          sequential=TRUE,
                          constraints = ~.,
                          monitor = NULL,
                          control = control.simulate.formula(),
                          verbose=FALSE, \dots)
\method{simulate}{ergm}(object, nsim=1, seed=NULL, coef=object$coef,
                       burnin=1000, interval=1000, 
                       statsonly=FALSE,
                       sequential=TRUE,
                       constraints = NULL,
                       monitor = NULL,
                       control = control.simulate.ergm(),
                       verbose=FALSE, \dots)
}
\description{
  \code{\link[stats]{simulate}} is used to draw from exponential family 
  random network models in their natural parameterizations.
  See \code{\link{ergm}} for more information on these models.
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
 \item{nsim}{Number of networks to be randomly drawn
    from the given distribution on the set of all networks,
    returned by the Metropolis-Hastings algorithm.}
 \item{seed}{Random number integer seed.
   The default is \code{sample(10000000, size=1)}.}
 \item{coef}{Vector of parameters for the model from which the sample is
   to be drawn.  If \code{object} is of class \code{ergm}, the default value
   is the vector of estimated coefficients.}
 \item{burnin}{The number of proposals before any MCMC sampling
    is done.  If NULL in simulate.ergm, the value of object$burnin is used.}
 \item{interval}{The number of proposals between sampled statistics.
    If NULL in simulate.ergm, the value of object$interval is used.
    The program prints a warning if too few proposals are being
    accepted (if the number of proposals between sampled observations
    ever equals an integral multiple of 100(1+the number of proposals
    accepted)).}
 \item{basis}{An optional \code{\link[network]{network}} object to start the MCMC
    algorithm from. This overrides the left-hand-side of the \code{formula}.
    If neither a left-hand-side nor a \code{basis} is present, an error results
    because the characteristics of the network (e.g., size and
    directedness) must be specified.}
  \item{constraints}{A one-sided formula specifying one or more constraints
    on the support of the distribution of the networks being
    simulated. See the documentation for a similar argument for
    \code{\link{ergm}} for more information. For
    \code{simulate.formula}, defaults to no constraints. For
    \code{simulate.ergm}, defaults to using the same constraints as
    those with which \code{object} was fitted.}
  \item{monitor}{A one-sided formula specifying one or more terms whose
    value is to be monitored. These terms are appeneded to the model, along
    with a coefficient of 0, so their statistics are returned.}
 \item{control}{A list of control parameters for algorithm
   tuning. Constructed using \code{\link{control.simulate.ergm}} or
   \code{\link{control.simulate.formula}}, which have different defaults.}
%\item{multiplicity}{The number of replications of proposal types
%  grouped together to form a single proposal.}
%\item{mixed}{Should a mixed parameterization be used for the
%  exponential family model with any degree terms in the model
%  using their mean-value parameters and the other terms using
%  the natural parameterization? If \code{FALSE} the natural
%  parameterization is used. If \code{FALSE} the coefficients
%  returned for the degree terms  are those that correspond to
%  the mean-values in the natural parameterization. The estimates
%  of the mean-values are the observed degree proportions.
%  }
%\item{fix}{Should the MCMC be run conditional on the non-degree 
%  terms? This gives the conditional estimates of these under
%  the natural parameterization, These are also the MLEs under
%  the mixed parameterization with the non-degree terms
%  represented by their mean-value parameters.
%  }
%\item{conditional}{Should a conditional parameterization be used for the
%  exponential family model with any degree terms in the model
%  using their mean-value parameters and the other terms using
%  the natural parameterization? If \code{FALSE} the natural
%  parameterization is used. If \code{TRUE} the coefficients
%  returned for the degree terms are their natural parameterization
%  estimates conditional on the non-degree term statistics.
%  }
 \item{statsonly}{If TRUE, return only the network statistics (not the 
   network(s) themselves)}
 \item{sequential}{If FALSE, each new simulation (of \code{nsim})
   is always started at the initial network.  If TRUE, the end of one simulation
   is used as the start of the next.
   Irrelevant when \code{nsim=1}.}
% \item{proposal}{The structure containing information about the
%  Metropolis-Hasting proposal used.}
 \item{verbose}{If this is \code{TRUE}, we will print out more information as 
  we run the program, including (currently) some goodness of fit 
  statistics.}
 \item{\dots}{further arguments passed to or used by methods.}
}

\details{
  A sample of networks is randomly drawn from the specified model.  The
  model is specified by the first argument of the function.  If
  the first argument is a \code{\link{formula}} then this defines 
  the model.  If the
  first argument is the output of a call to \code{\link{ergm}} then the model
  used for that call is the one fit - and unless \code{coef} is specified,
  the sample is from the MLE of the parameters.  
  If neither of those are given as the
  first argument then a Bernoulli network is generated with the probability
  of ties defined by \code{prob} or \code{coef}.

  Note that the first network is sampled after \code{burnin} + \code{interval}
  steps, and any subsequent networks are sampled each
  \code{interval} steps after the first.

  More information can be found by looking at the documentation of
  \code{\link{ergm}}.
}

\value{If \code{nsim==1}, \code{simulate.ergm} returns an object of 
class \code{network}.  If \code{nsim>1}, it returns an object of class
\code{network.list} that is a list with the following elements:
  \item{formula}{The \code{\link{formula}} used to generate the sample.}
  \item{networks}{A list of the generated networks.}
  \item{stats}{The \eqn{n\times p}  matrix of network change statistics, 
   where \eqn{n} is the sample size and \eqn{p} is the number
   of network change statistics specified in the model.}
}
\seealso{ergm, network, print.network}
\examples{
#
# Let's draw from a Bernoulli model with 16 nodes
# and density 0.5 (i.e., coef = c(0,0))
#
g.sim <- simulate(network(16) ~ edges + mutual, coef=c(0, 0))
#
# What are the statistics like?
#
summary(g.sim ~ edges + mutual)
#
# Now simulate a network with higher mutuality
#
g.sim <- simulate(network(16) ~ edges + mutual, coef=c(0,2))
#
# How do the statistics look?
#
summary(g.sim ~ edges + mutual)
#
# Let's draw from a Bernoulli model with 16 nodes
# and tie probability 0.1
#
g.use <- network(16,density=0.1,directed=FALSE)
#
# Starting from this network let's draw 5 realizations
# of a edges and 2-star network
#
g.sim <- simulate(~edges+kstar(2), nsim=5, coef=c(-1.8,0.03),
               basis=g.use, control=control.simulate(
                 MCMC.burnin=100000,
                 MCMC.interval=1000))
g.sim
summary(g.sim)
#
# attach the Florentine Marriage data
#
data(florentine)
#
# fit an edges and 2-star model using the ergm function
#
gest <- ergm(flomarriage ~ edges + kstar(2))
summary(gest)
#
# Draw from the fitted model (satatistics only), and observe the number
# of triangles as well.
#
g.sim <- simulate(gest, nsim=100, 
            monitor=~triangles, statsonly=TRUE,
            control=control.simulate(MCMC.burnin=1000, MCMC.interval=1000))
head(g.sim)
}
\keyword{models}
