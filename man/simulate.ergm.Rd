% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.ergm.R
\name{simulate.ergm}
\alias{simulate.ergm}
\alias{simulate.formula_lhs_network}
\alias{simulate.formula.ergm}
\alias{simulate_formula}
\alias{.simulate_formula.network}
\alias{simulate_formula.network}
\alias{simulate_formula.ergm_state}
\alias{simulate.ergm_model}
\alias{simulate.ergm_state_full}
\alias{simulate.ergm_state}
\title{Draw from the distribution of an Exponential Family Random Graph Model}
\usage{
\method{simulate}{formula_lhs_network}(object, nsim = 1, seed = NULL, ...)

simulate_formula(object, ..., basis = eval_lhs.formula(object))

\method{simulate_formula}{network}(
  object,
  nsim = 1,
  seed = NULL,
  coef,
  response = NULL,
  reference = ~Bernoulli,
  constraints = ~.,
  observational = FALSE,
  monitor = NULL,
  statsonly = FALSE,
  esteq = FALSE,
  output = c("network", "stats", "edgelist", "ergm_state"),
  simplify = TRUE,
  sequential = TRUE,
  control = control.simulate.formula(),
  verbose = FALSE,
  ...,
  basis = ergm.getnetwork(object),
  do.sim = NULL,
  return.args = NULL
)

\method{simulate_formula}{ergm_state}(
  object,
  nsim = 1,
  seed = NULL,
  coef,
  response = NULL,
  reference = ~Bernoulli,
  constraints = ~.,
  observational = FALSE,
  monitor = NULL,
  statsonly = FALSE,
  esteq = FALSE,
  output = c("network", "stats", "edgelist", "ergm_state"),
  simplify = TRUE,
  sequential = TRUE,
  control = control.simulate.formula(),
  verbose = FALSE,
  ...,
  basis = ergm.getnetwork(object),
  do.sim = NULL,
  return.args = NULL
)

\method{simulate}{ergm_model}(
  object,
  nsim = 1,
  seed = NULL,
  coef,
  reference = if (is(constraints, "ergm_proposal")) NULL else trim_env(~Bernoulli),
  constraints = trim_env(~.),
  observational = FALSE,
  monitor = NULL,
  basis = NULL,
  esteq = FALSE,
  output = c("network", "stats", "edgelist", "ergm_state"),
  simplify = TRUE,
  sequential = TRUE,
  control = control.simulate.formula(),
  verbose = FALSE,
  ...,
  do.sim = NULL,
  return.args = NULL
)

\method{simulate}{ergm_state_full}(
  object,
  nsim = 1,
  seed = NULL,
  coef,
  esteq = FALSE,
  output = c("network", "stats", "edgelist", "ergm_state"),
  simplify = TRUE,
  sequential = TRUE,
  control = control.simulate.formula(),
  verbose = FALSE,
  ...,
  return.args = NULL
)

\method{simulate}{ergm}(
  object,
  nsim = 1,
  seed = NULL,
  coef = coefficients(object),
  response = object$network \%ergmlhs\% "response",
  reference = object$reference,
  constraints = list(object$constraints, object$obs.constraints),
  observational = FALSE,
  monitor = NULL,
  basis = if (observational) object$network else NVL(object$newnetwork, object$network),
  statsonly = FALSE,
  esteq = FALSE,
  output = c("network", "stats", "edgelist", "ergm_state"),
  simplify = TRUE,
  sequential = TRUE,
  control = control.simulate.ergm(),
  verbose = FALSE,
  ...,
  return.args = NULL
)
}
\arguments{
\item{object}{Either a \code{\link{formula}} or an
\code{\link{ergm}} object.  The \code{\link{formula}} should be of the form
\code{y ~ <model terms>}, where \code{y} is a network object or a matrix
that can be coerced to a \code{\link[network:network]{network}} object.  For the
details on the possible \code{<model terms>}, see \code{\link{ergmTerm}}.
To create a \code{\link[network:network]{network}} object in , use the
\code{network()} function, then add nodal attributes to it using the
\code{\%v\%} operator if necessary.}

\item{nsim}{Number of networks to be randomly drawn from the given
distribution on the set of all networks, returned by the Metropolis-Hastings
algorithm.}

\item{seed}{Seed value (integer) for the random number generator.  See
\code{\link[=set.seed]{set.seed()}}.}

\item{\dots}{Further arguments passed to or used by methods.}

\item{basis}{a value (usually a \code{\link[network:network]{network}}) to override the LHS of the formula.}

\item{coef}{Vector of parameter values for the model from which the
sample is to be drawn.  If \code{object} is of class \code{ergm},
the default value is the vector of estimated coefficients. Can be
set to \code{NULL} to bypass, but only if \code{return.args} below is used. \matchnames{coefficient}}

\item{response}{Either a character string, a formula, or \code{NULL} (the default), to specify the response attributes and whether the ERGM is binary or valued. Interpreted as follows: \describe{

\item{\code{NULL}}{Model simple presence or absence, via a binary ERGM.}

\item{character string}{The name of the edge attribute whose value is to be modeled. Type of ERGM will be determined by whether the attribute is \code{\link{logical}} (\code{TRUE}/\code{FALSE}) for binary or \code{\link{numeric}} for valued.}

\item{a formula}{must be of the form \code{NAME~EXPR|TYPE} (with \code{|} being literal). \code{EXPR} is evaluated in the formula's environment with the network's edge attributes accessible as variables. The optional \code{NAME} specifies the name of the edge attribute into which the results should be stored, with the default being a concise version of \code{EXPR}. Normally, the type of ERGM is determined by whether the result of evaluating \code{EXPR} is logical or numeric, but the optional \code{TYPE} can be used to override by specifying a scalar of the type involved (e.g., \code{TRUE} for binary and \code{1} for valued).}

}}

\item{reference}{A one-sided formula specifying
the reference measure (\eqn{h(y)}) to be used.
See help for \link[=ergm-references]{ERGM reference measures} implemented in the
\strong{\link[=ergm-package]{ergm}} package.}

\item{constraints}{A formula specifying one or more constraints
on the support of the distribution of the networks being modeled. Multiple constraints
may be given, separated by \dQuote{+} and \dQuote{-} operators. See
\code{\link{ergmConstraint}} for the detailed explanation of
their semantics and also for an indexed list of the constraints visible to the \pkg{ergm} package.

The default is to have no constraints except those provided through
the \code{\link{ergmlhs}} API.

Together with the model terms in the formula and the reference measure, the constraints
define the distribution of networks being modeled.

It is also possible to specify a proposal function directly either
by passing a string with the function's name (in which case,
arguments to the proposal should be specified through the
\code{MCMC.prop.args} argument to the relevant control function, or
by giving it on the LHS of the hints formula to \code{MCMC.prop}
argument to the control function. This will override
the one chosen automatically.

Note that not all possible combinations of constraints and reference
measures are supported. However, for relatively simple constraints
(i.e., those that simply permit or forbid specific dyads or sets of
dyads from changing), arbitrary combinations should be possible.}

\item{observational}{Inherit observational constraints rather than model
constraints.}

\item{monitor}{A one-sided formula specifying one or more terms
whose value is to be monitored. These terms are appended to the
model, along with a coefficient of 0, so their statistics are
returned. An \code{\link{ergm_model}} objectcan be passed as well.}

\item{statsonly}{Logical: If TRUE, return only the network statistics, not
the network(s) themselves. Deprecated in favor of \verb{output=}.}

\item{esteq}{Logical: If TRUE, compute the sample estimating equations of an
ERGM: if the model is non-curved, all non-offset statistics are returned
either way, but if the model is curved, the score estimating function values
(3.1) by Hunter and Handcock (2006) are returned instead.}

\item{output}{Normally character, one of \code{"network"} (default),
\code{"stats"}, \code{"edgelist"}, or \code{"ergm_state"}: determines the output
format. Partial matching is performed.

Alternatively, a function with prototype
\verb{function(ergm_state, chain, iter, ...)} that is
called for each returned network, and its return value, rather
than the network itself, is stored. This can be used to, for
example, store the simulated networks to disk without storing
them in memory or compute network statistics not implemented
using the ERGM API, without having to store the networks
themselves.}

\item{simplify}{Logical: If \code{TRUE} the output is "simplified":
sampled networks are returned in a single list, statistics from
multiple parallel chains are stacked, etc.. This makes it
consistent with behavior prior to \code{ergm} 3.10.}

\item{sequential}{Logical: If FALSE, each of the \code{nsim} simulated
Markov chains begins at the initial network.  If TRUE, the end of one
simulation is used as the start of the next.  Irrelevant when \code{nsim=1}.}

\item{control}{A list of control parameters for algorithm tuning,
typically constructed with \code{\link[=control.simulate.ergm]{control.simulate.ergm()}} or \code{\link[=control.simulate.formula]{control.simulate.formula()}}, which have different
defaults. Their documentation gives the the list of recognized
control parameters and their meaning. The more generic utility
\code{\link[=snctrl]{snctrl()}} (StatNet ConTRoL) also provides argument completion
for the available control functions and limited argument name
checking.}

\item{verbose}{A logical or an integer to control the amount of
progress and diagnostic information to be printed. \code{FALSE}/\code{0}
produces minimal output, with higher values producing more
detail. Note that very high values (5+) may significantly slow
down processing.}

\item{do.sim}{Logical; a deprecated interface superseded by \code{return.args},
that saves the inputs to the next level of the function.}

\item{return.args}{Character; if not \code{NULL}, the \code{simulate} method
for that particular class will, instead of proceeding for
simulation, instead return its arguments as a list that can be
passed as a second argument to \code{\link[=do.call]{do.call()}} or a lower-level
function such as \code{\link[=ergm_MCMC_sample]{ergm_MCMC_sample()}}. This can be useful if, for
example, one wants to run several simulations with varying
coefficients and does not want to reinitialize the model and the
proposal every time. Valid inputs at this time are \code{"formula"},
"ergm_model", and one of the \code{"ergm_state"} classes, for the three
respective stopping points.}
}
\value{
If \code{output=="stats"} an \code{\link[coda:mcmc]{mcmc}} object containing the
simulated network statistics. If \code{control$parallel>0}, an
\code{\link[coda:mcmc.list]{mcmc.list}} object. If \code{simplify=TRUE} (the default), these
would then be "stacked" and converted to a standard \code{\link{matrix}}. A
logical vector indicating whether or not the term had come from
the \verb{monitor=} formula is stored in \code{\link[=attr]{attr()}}-style attribute
\code{"monitored"}.

Otherwise, a representation of the simulated network is returned,
in the form specified by \code{output}. In addition to a network
representation or a list thereof, they have the following
\code{\link[=attr]{attr()}}-style attributes: \describe{

\item{\code{formula}}{The \code{\link{formula}} used to generate the
sample.}

\item{\code{stats}}{An \code{\link[coda:mcmc]{mcmc}} or \code{\link[coda:mcmc.list]{mcmc.list}} object as above.}

\item{\code{control}}{Control parameters used to generate the sample.}

\item{\code{constraints}}{Constraints used to generate the sample.}

\item{\code{reference}}{The reference measure for the sample.}

\item{\code{monitor}}{The monitoring formula.}

\item{\code{response}}{The edge attribute used as a response.}

}

The following are the permitted network formats: \describe{

\item{\code{"network"}}{A \code{\link[network:network]{network}} object.}

\item{\code{"edgelist"}}{An \code{\link[network:as.edgelist]{edgelist}} representation of the network,
or a list thereof, depending on \code{nsim}.}

\item{\code{"ergm_state"}}{\code{\link{ergm_state}}, a semi-internal representation of
a network consisting of a \code{\link[network:network]{network}} object emptied of edges, with
an attached edgelist matrix.}

}

If \code{nsim > 1}, these objects are returned in a list, with
class \code{\link{network.list}}.

If \code{simplify==FALSE}, the networks are returned as a nested list,
with outer list being the parallel chain (including 1 for no
parallelism) and inner list being the samples within that chains
(including 1, if one network per chain). It is assigned an
additional class, \code{network.list.list}. If \code{TRUE}, they are
concatenated, and if a total of one network had been simulated, the
network itself will be returned.
}
\description{
\code{\link[stats]{simulate}} is used to draw from exponential
family random network models.  See \code{\link[=ergm]{ergm()}} for more
information on these models.

The method for \code{\link{ergm}} objects inherits the model,
the coefficients, the response attribute, the reference, the
constraints, and most simulation parameters from the model fit,
unless overridden by passing them explicitly. Unless overridden,
the simulation is initialized with either a random draw from near
the fitted model saved by \code{\link[=ergm]{ergm()}} or, if unavailable, the
network to which the ERGM was fit.
}
\details{
A sample of networks is randomly drawn from the specified model.  The model
is specified by the first argument of the function.  If the first argument
is a \code{\link{formula}} then this defines the model.  If the first
argument is the output of a call to \code{\link[=ergm]{ergm()}} then the model used
for that call is the one fit -- and unless \code{coef} is specified, the
sample is from the MLE of the parameters.  If neither of those are given as
the first argument then a Bernoulli network is generated with the
probability of ties defined by \code{prob} or \code{coef}.

Note that the first network is sampled after \code{burnin} steps,
and any subsequent networks are sampled each \code{interval} steps
after the first.

More information can be found by looking at the documentation of
\code{\link[=ergm]{ergm()}}.
}
\section{Functions}{
\itemize{
\item \code{simulate(ergm_state_full)}: a low-level function to simulate from an \code{\link{ergm_state}} object.

}}
\note{
The actual \code{\link[network:network]{network}} method for \code{\link[=simulate_formula]{simulate_formula()}} is
actually called \code{.simulate_formula.network()} and is also
exported as an object. This allows it to be overridden by
extension packages, such as \code{tergm}, but also accessed directly
when needed.

\code{\link[=simulate.ergm_model]{simulate.ergm_model()}} is a lower-level interface, providing
a \code{\link[=simulate]{simulate()}} method for the \code{\link{ergm_model}} class. The \code{basis}
argument is required; \code{monitor}, if passed, must be an
\code{\link{ergm_model}} as well; and \code{constraints} can be an
\code{\link{ergm_proposal}} object instead.
}
\examples{
\dontshow{
options(ergm.eval.loglik=FALSE)
}
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
# Starting from this network let's draw 3 realizations
# of a edges and 2-star network
#
g.sim <- simulate(~edges+kstar(2), nsim=3, coef=c(-1.8,0.03),
               basis=g.use, control=control.simulate(
                 MCMC.burnin=1000,
                 MCMC.interval=100))
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
# Draw from the fitted model (statistics only), and observe the number
# of triangles as well.
#
g.sim <- simulate(gest, nsim=10, 
            monitor=~triangles, output="stats",
            control=control.simulate.ergm(MCMC.burnin=1000, MCMC.interval=100))
g.sim

# Custom output: store the edgecount (computed in R), iteration index, and chain index.
output.f <- function(x, iter, chain, ...){
  list(nedges = network.edgecount(as.network(x)),
       chain = chain, iter = iter)
}
g.sim <- simulate(gest, nsim=3,
            output=output.f, simplify=FALSE,
            control=control.simulate.ergm(MCMC.burnin=1000, MCMC.interval=100))
unclass(g.sim)
}
\seealso{
\code{\link[=ergm]{ergm()}}, \code{\link[network:network]{network}},
\code{\link[=ergm_MCMC_sample]{ergm_MCMC_sample()}} for a demonstration of \verb{return.args=}.
}
\keyword{models}
