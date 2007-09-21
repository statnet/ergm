\name{ergm}
\alias{ergm}
\title{Exponential Family Random Graph Models}
\description{
  \code{\link{ergm}} is used to fit linear exponential random network models, in which
  the probability of a given network, \eqn{y}, on a set of nodes is 
  \eqn{\exp(\theta{\cdot}g(y))/c(\theta)}, where 
  \eqn{g(y)} is a vector of network statistics,
  \eqn{\theta} is a parameter vector of the same length and \eqn{c(\theta)} is the
  normalizing constant for the distribution.
  \code{\link{ergm}} can return either a maximum pseudo-likelihood
  estimate or an approximate maximum likelihood estimator based on a Monte Carlo scheme.
}
\usage{
ergm(formula, theta0="MPLE",
     MPLEonly=FALSE, MLestimate=!MPLEonly, seed=NULL,
     burnin=10000, MCMCsamplesize=10000, interval=100, maxit=3,
     proposaltype="randomtoggle",
     algorithm.control=list(),
     verbose=FALSE, \dots)
}
\arguments{
 \item{formula}{formula; an \R \code{\link{formula}} object, of the form
 \code{y ~ <model terms>},
  where \code{y} is a network object or a matrix that can be
  coerced to a network object.  For the details on the possible
  \code{<model terms>}, see \code{\link{ergm-terms}}.  To create a
  network object in \R, use the \code{network()} function,
  then add nodal attributes to it using the \code{\%v\%}
  operator if necessary.}
 \item{theta0}{vector; the parameter value used to generate
  the MCMC sample and as a starting value for the estimation. By default the MPLE is used (\code{theta0="MPLE"}).}
 \item{MPLEonly}{logical; \code{TRUE} if the maximum
  pseudo-likelihood estimate is to be computed and returned.
  Note that \code{MPLEonly=TRUE} will render moot most other parameters in this list.}
 \item{MLestimate}{logical; \code{TRUE} if only the Monte Carlo maximum likelihood estimate
  is to be computed and returned.}
 \item{burnin}{count; the number of proposals before any MCMC sampling
  is done. It typically is set to a fairly large number.}
%  Currently, there is no support for any check of the Markov
%  chain mixing, so \code{burnin} should be set to a fairly large number.}
 \item{MCMCsamplesize}{count; the number of network statistics,
  randomly drawn from a given distribution on the set of all networks,
  returned by the Metropolis-Hastings algorithm.}
 \item{interval}{count; the number of proposals between sampled statistics.
% The program prints a warning if too few proposals are being
% accepted in any interval before each sample.
% (if the number of proposals between sampled observations
% ever equals an integral multiple of 100(1+the number of proposals accepted)).
  }
 \item{maxit}{count; the number of times the parameter for the MCMC
  should be updated by maximizing the MCMC likelihood. At each
  step the parameter is changed to the values that maximizes
  the MCMC likelihood based on the current sample.}
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
 \item{algorithm.control}{A list of control parameters. See 
  \bold{Algorithm Tuning}.}
 \item{seed}{integer; random number integer seed.
   The default is \code{sample(10000000, size=1)}.}
 \item{verbose}{logical; if this is \code{TRUE}, the program will print out additional information, including goodness of fit 
  statistics.}
 \item{\dots}{Additional arguments, to be passed to lower-level functions
  in the future.}
}

\section{Model Terms}{
  The \code{\link{ergm}} function allows the user to explore a large number
  of potential models for their network data.  The
  terms currently supported by the program,
  and a brief description of each is given in the documentation
  \code{\link{ergm-terms}}.
  In the \code{\link{formula}} for the model, the model terms are various function-like
  calls, some of which require arguments, separated by \code{+} signs.
  See \code{\link{ergm-terms}} for details.
}

\value{\code{\link{ergm}} returns an object of class \code{\link{ergm}} that is a list
  consisting of the following elements:
 \item{coef}{The Monte Carlo maximum likelihood estimate
  of \eqn{\theta}, the vector of coefficients for the model
  parameters.}
 \item{sample}{The \eqn{n\times p} matrix of network statistics, 
  where \eqn{n} is the
  sample size and \eqn{p} is the number of network statistics specified in the
  model, that is used in the maximum likelihood estimation routine.}
 \item{iterations}{The number of Newton-Raphson iterations required
  before convergence.}
 \item{MCMCtheta}{The value of \eqn{\theta} used to produce the Markov chain
  Monte Carlo sample.  As long as the Markov chain mixes sufficiently
  well, \code{sample} is roughly a random sample from the distribution
  of network statistics specified by the model with the parameter equal
  to \code{MCMCtheta}.  If \code{MPLEonly=TRUE} then 
  \code{MCMCtheta} equals the MPLE.}
 \item{loglikelihood}{The approximate log-likelihood for the MLE.
  The value is only approximate because it is estimated based 
  on the MCMC random sample.}
 \item{gradient}{The value of the gradient vector of the approximated
  loglikelihood function, evaluated at the maximizer.  This vector
  should be very close to zero.}
 \item{hessian}{The Hessian matrix of the approximated loglikelihood
  function, evaluated at the maximizer.  This matrix may be inverted to
  give an approximate covariance matrix for the MLE.}
 \item{samplesize}{The size of the MCMC sample}
 \item{formula}{The original \code{\link{formula}} entered into the \code{\link{ergm}} function.}
 \item{statsmatrix}{If the option \code{returnMCMCstats=TRUE}, this is 
  the the matrix of change statistics from the MCMC run.}
 \item{newnetwork}{The network generated at the end of the MCMC sampling.}

  See the method \code{\link{print.ergm}} for details on how
  an \code{\link{ergm}} object is printed.  Note that the
  method \code{\link{summary.ergm}} returns a summary of the
  relevant parts of the \code{\link{ergm}} object in concise summary
  format.
}

\section{Notes on model specification}{
Although each of the statistics in a given model is a summary
statistic for the entire network, it is rarely necessary to
calculate statistics for an entire network
in a proposed Metropolis-Hastings step.
%  The reason for this is that we may assume that each statistic equals
% zero for the observed network (translating each statistic by
% a constant does not change the probability model), and then
% it is only necessary to calculate the changes in the various
%statistics for a proposed Metropolis-Hastings step.  
Thus, for example, if the triangle term is included in the model,
a census of all triangles in the observed network is never
taken; instead, only the change in the number of triangles
is recorded for each edge toggle.

In the implementation of \code{\link{ergm}}, the model is initialized
in \R, then all the model information is passed to a C program
that generates the sample of network statistics using MCMC.
This sample is then returned to \R, which implements a
simple Newton-Raphson algorithm to approximate the MLE.
An alternative style of maximum likelihood estimation is to use a stochastic
approximation algorithm. This can be chosen with the 
\code{algorithm.control(style="Robbins-Monro")} option.
%There seems little reason to implement the Newton-Raphson
%algorithm in C since it moves very quickly in \R and
%it is easier to alter the code if it is written in \R.

The default mechanism for proposing new networks for the MCMC sample 
space is the Metropolis-Hastings algorithm, which simply chooses a dyad 
at random and proposes to toggle that edge; each possible dyad is 
equally likely.  The \code{proposaltype} 
option allow many more complex proposals to be specified.
We have developed and implemented a wide range of algorithms. These are described
in the documentation for \code{proposaltype}.
%, including a double
%Metropolis-Hastings algorithm (two random toggles); an algorithm
%that reallocates all the edges of one node, with replacement;
%and an algorithm that toggles all dyads based on a single node
%(i.e., gives the complement edges for all dyads rooted at
%a randomly chosen node).
For example, we have included proposal functions that condition on
maintaining the absolute degree distribution for the observed
network.  Each proposal network will have exactly the same number
of nodes with each degree as does the original network; this
means that if the proposal network removes an edge between a
node of degree 3 and a node of degree 5, it must also add
an edge between a node of degree 2 and a node of degree 4.
Note that one or both of the latter nodes may be the same as
the former nodes.
%See the \code{proposaltype} argument for details.

The package is designed so that the user can add additional 
proposal types. 
}
  

\section{Conditioning within models}{
There are many times when one may wish to condition on the
number of inedges or outedges possessed by a node, either as a
consequence of some intrinsic property of that node (e.g., to control for 
activity or popularity processes), to account
for known outliers of some kind, and thus we wish to limit its indegree, an
intrinsic property of the sampling scheme whence came our data (e.g.,
the survey asked everyone to name only three friends total) or as a
function of the attributes of the nodes to which a node has edges
(e.g., we specify that nodes designated ``male'' have a maximum number
of outdegrees to nodes designated ``female'').  To accomplish this we
use the \code{\link{ergm}} parameter \code{boundDeg}.

Let's consider the simple cases first. 
Suppose you want to condition on the total number of degrees
regardless of attributes.  That is, if you had a survey that asked
respondents to name three alters and no more, then you might want to
limit your maximal outdegree to three without regard to any of the
alters' attributes. The argument is then:  

   \code{boundDeg=list(maxout=3)}

Similar calls are used to restrict the number of indegrees
(\code{maxin}), the minimum number of outdegrees
(\code{minout}), and the minimum number of indegrees
(\code{minin}).

You can also set ego specific limits. For example:

   \code{boundDeg=list(maxout=rep(c(3,4),c(36,35)))}

limits the first 36 to 3 and the other 35 to 4 outdegrees.

Multiple restrictions can be combined. \code{boundDeg} is very flexible.
In general, the \code{boundDeg} term can contain up to six arguments:

  \preformatted{
 boundDeg <- list(attribs=attribs,
	          maxout=maxout,
                  maxin=maxin,
                  minout=minout,
                  minin=minin,
                  condAllDegExact=FALSE)
  }

Omitted arguments are unrestricted, and arguments of length 1
are replicated out to all nodes (as above). If an individual
entry in \code{maxout},..., \code{minin} is \code{NA} then
no restriction of that kind is applied to that actor.

In general, \code{attribs} is a matrix of the attributes on
which we are conditioning.  The dimensions of \code{attribs}
are \code{n_nodes} rows by \code{attrcount} columns, where
\code{attrcount} is the number of distinct attribute values
on which we want to condition (i.e., a separate column is
required for ``male'' and ``female'' if we want to condition on
the number of ties to both ``male'' and ``female'' partners).
The value of \code{attribs[n, i]}, therefore, is \code{TRUE}
if node \code{n} has attribute value i, and \code{FALSE} otherwise.
(Note that, since each column represents only a single value
of a single attribute, the values of this matrix are all
Boolean (\code{TRUE} or \code{FALSE}).)  It is important to
note that \code{attribs} is a matrix of nodal attributes,
not alter attributes.

So, for instance, if we wanted to construct an \code{attribs} matrix
with two columns, one each for male and female attribute
values (we are conditioning on these values of the attribute
``sex''), and the attribute sex is represented in ads.sex as
an \code{n_node}-long vector of 0s and 1s (men and women),
then our code would look as follows:

\preformatted{
 # male column: bit vector, TRUE for males
 attrsex1 <- (ads.sex == 0) 
 # female column: bit vector, TRUE for females
 attrsex2 <- (ads.sex == 1)
 # now create attribs matrix
 attribs <- matrix(ncol=2,nrow=71, data=c(attrsex1,attrsex2))
}

\code{maxout} is a matrix of alter attributes, with the same
dimensions as the \code{attribs} matrix.  \code{maxout} is \code{n_nodes}
rows by \code{attrcount} columns.  The value of \code{maxout[n,i]},
therefore, is the maximum number of outdegrees permitted
from node \code{n} to nodes with the attribute \code{i} (where a \code{NA}
means there is no maximum).

For example:  if we wanted to create a \code{maxout} matrix to work
with our \code{attribs} matrix above, with a maximum from every
node of five outedges to males and five outedges to females,
our code would look like this:

\preformatted{
      # every node has maximum of 5 outdegrees to male alters
      maxoutsex1 <- c(rep(5,71))
      # every node has maximum of 5 outdegrees to female alters
      maxoutsex2 <- c(rep(5,71))
      # now create maxout matrix
      maxout <- cbind(maxoutsex1,maxoutsex2)
}

The \code{maxin}, \code{minout}, and \code{minin} matrices
are constructed exactly like the \code{maxout} matrix,
except for the maximum allowed indegree, the minimum allowed
outdegree, and the minimum allowed indegree, respectively.
Note that in an undirected network, we only look at the outdegree
matrices; \code{maxin} and \code{minin} will both be ignored
in this case.

% First, you would assign the ``ego'' attribute to all your nodes -
% because every node is an ego, every node gets a \code{TRUE} (or 1) in its
% column in the \code{attribs} matrix.  

% \code{attribs <- matrix(ncol=1,nrow=71,data=c(rep(1,71)))}

% Then you would assign a maximum number of alter out degrees allowed
% to each ego in the same column as in the \code{attribs} matrix.  In this
% case, all egos have a maximum outdegree of 3.

% \code{maxout <- matrix(ncol=1,nrow=71,data=c(rep(3,71)))}

% Then, we can create dummy \code{maxin}, \code{minout},
% \code{minin} matrices (in the future these will be empty
% by default) and create the \code{boundDeg} argument to pass
% into the \code{\link{ergm}} call.

% To create ``dummy'' \code{maxin}, \code{minout}, and
% \code{minin} matrices to correspond to the \code{attribs} and
% \code{maxout} matrices created above, that impose no limits
% on the amount of edges allowed, our code would look like this:

% \code{maxin <- matrix(NA,ncol=2,nrow=71)
%       minout <- matrix(NA,ncol=2,nrow=71)
%	minin <- matrix(NA,ncol=2,nrow=71)}

%  The argument we pass into the \code{\link{ergm}} call, then would look like this:

%  \code{boundDeg = list( attribs=attribs,
%        	         maxout=maxout,
%        	         maxin=maxin,
%        	         minout=minout,
%        	         minin=minin)}

%  This argument passed into \code{\link{ergm}} will now restrict your proposed 
%  networks to those which meet your desired degree conditions.
  

%  \code{boundDeg = list( attribs=attribs,
%        	         maxout=maxout,
%        	         maxin=matrix(NA,ncol=1,nrow=71),
%        	         minout=matrix(NA,ncol=1,nrow=71),
%        	         minin=matrix(NA,ncol=1,nrow=71))}

  Finally, \code{boundDeg} takes an optional argument,
  \code{condAllDegExact}.  If \code{condAllDegExact} is
  \code{TRUE}, then all other inputs to this function will be
  ignored, and we will condition on the initial degree of all
  nodes in the observed network - no node will change degree in
  any toggle.  To achieve this, the code automatically creates
  the five matrices as if only the ego column were specified.
  It then assigns values to the ego column, for each \eqn{n} in
  \eqn{1\dots\code{n_nodes}}.

  \preformatted{
attribs[n][0] = 1 # just the ego values
maxout[n][0] = minout[n][0] = observed outdegree of n in network
maxin[n][0] = minin[n][0] = observed indegree of n in network
  }

% Note that if \code{condAllDegExact} is \code{TRUE}, then
% \code{toggletwo} must be \code{TRUE} as well or else there will never
% be a valid toggle.
}

\section{Algorithm Tuning}{
  
  There are a large number of parameters that can be modified to define complex
  models, estimators and computational aspects. These are specified via the 
  \code{algorithm.control} argument.
  The \code{algorithm.control} argument is a list that can supply any of the
  following components:
  \describe{
    \item{\code{nr.maxit}}{count; The maximum number of iterations in the
      Newton-Raphson optimization.  Defaults to \code{100}.
      \code{maxit} gives the total number of likelihood
      function evaluations.}
    \item{\code{calc.mcmc.se}}{logical; should the contribution to the 
      standard errors of the estimator incurred by the MCMC sampling
      be computed. Default is \code{TRUE}.}
    \item{\code{hessian}}{logical; Should the Hessian matrix
      of the likelihood function be computed. 
      Default is \code{TRUE}.}
    \item{\code{compress}}{logical; Should the matrix of sample statistics
      returned be compressed to the set of unique statistics with a 
      column of frequencies post-pended.  This also uses a compression
      algorithm in the computation of the maximum psuedo-likelihood
      estimate that will dramatically speed it for large networks.
      Default is \code{FALSE}.}
    \item{\code{maxNumDyadTypes}}{count; The maximum number of unique
      pseudolikelihood change statistics to be allowed if \code{compress=TRUE}.
      It is only relevant in that case.
      Default is \code{10000}.}
    \item{\code{MPLEsamplesize}}{count; the sample size to use for endogenous
      sampling in the pseudolikelihood computation.
      Default is \code{50000}.}
    \item{\code{trace}}{non-negative integer; If positive,
      tracing information on the
      progress of the optimization is produced. Higher values may
      produce more tracing information: for method \code{"L-BFGS-B"}
      there are six levels of tracing.  (To understand exactly what
      these do see the source code for \code{\link[stats]{optim}}: higher levels 
      give more detail.)}
    \item{\code{boundDeg}}{list; comprising up to five matrices and 
      one or more arguments, which are used to constrain the number of in 
      or out degrees that a node can have. By default no constraining 
      is applied.  This is explained more fully under the 
      \code{Conditioning within models} section. Defaults to no conditioning.}
    \item{\code{drop}}{logical; Should the degenerate terms in the model be
      dropped from the fit?
      If statistics occur on the extreme of their range they
      correspond to infinite parameter estimates.
      Default is \code{TRUE}.}
    \item{\code{proposalpackage}}{character; The name of the package the 
      \code{proposaltype} is from.
      Defaults to \code{"ergm"}.}
    \item{\code{mcmc.precision}}{vector; upper bounds on the precision of the 
      standard errors induced by the MCMC algorithm.
      Defaults to \code{0.05}.}
    \item{\code{metric}}{character; The name of the optimization metric
      to use. Defaults to \code{"Likelihood"}.}
    \item{\code{method}}{character; The name of the optimization method
      to use. See \code{\link[stats]{optim}} for the options. The default method
      \code{"BFGS"} is a quasi-Newton method (also known as a variable
      metric algorithm). It is attributed to
      Broyden, Fletcher, Goldfarb and Shanno. This uses function values
      and gradients to build up a picture of the surface to be optimized.}
    \item{\code{style}}{character; The style of maximum 
      likelihood estimation to use. The default is optimization of an
      MCMC estimate of the log-likelihood. An alternative is to use 
      a form of stochastic approximation (\code{"Robbins-Monro"}).
      The direct use of the likelihood function has many theoretical
      advantages over stochastic approximation, but the choice will
      depend on the model and data being fit. See Handcock (2000) and 
      Hunter and Handcock (2006) for details.}
    \item{\code{trustregion}}{numeric; The maximum amount the algorithm will
      allow the approximated likelihood to be increased at a given iteration.
      Defaults to 20.
      See Boer, Huisman, Snijders, and Zeggelink (2003) for details.}
    \item{\code{phase1_n}}{count; The number of MCMC samples to draw
      in Phase 1 of the stochastic approximation algorithm.
      Defaults to 7 plus 3 times the number of terms in the model.
      See Boer, Huisman, Snijders, and Zeggelink (2003) for details.}
    \item{\code{nsubphases}}{count; The number of sub-phases 
      in Phase 2 of the stochastic approximation algorithm.
      Defaults to \code{maxit}.
      See Boer, Huisman, Snijders, and Zeggelink (2003) for details.}
    \item{\code{initial_gain}}{numeric; The initial gain to
      Phase 2 of the stochastic approximation algorithm.
      Defaults to 0.1.
      See Boer, Huisman, Snijders, and Zeggelink (2003) for details.}
    \item{\code{niterations}}{count; The number of MCMC samples to draw
      in Phase 2 of the stochastic approximation algorithm.
      Defaults to 7 plus the number of terms in the model.
      See Boer, Huisman, Snijders, and Zeggelink (2003) for details.}
    \item{\code{phase3_n}}{count; The sample size for the MCMC sample
      in Phase 3 of the stochastic approximation algorithm.
      Defaults to 1000.
      See Boer, Huisman, Snijders, and Zeggelink (2003) for details.}
    \item{\code{returnMCMCstats}}{logical; If this is \code{TRUE} (the
      default) the matrix of change 
      statistics from the MCMC run is returned as component \code{sample}.
      This matrix is actually an object of class \code{mcmc} and can be 
      used directly in the \code{CODA} package to assess MCMC
      convergence.}
  }
}

 
\section{Dealing with degeneracy}{
In order to begin the process of estimating network coefficients,
we need starting values - guesses at the true values of the
network statistic coefficients.  The default is to begin the MCMC
estimation process at the deterministic MPLE values. These
values are often taken as good-enough final answers by many other
applications. However recent work has indicated that they are sub-optimal and
can be dramatically bad.

In using the MPLE values, MCMC MLE often runs into
problems caused by the inherent instability of the natural parameter
space of the models (Handcock 2000, 2002, 2003).
If the initial values for the parameter coefficients
are off by a very small amount in the wrong direction, the
result is often a sample of networks that are degenerate -
that is, networks that are entirely full or entirely empty,
or that are otherwise less-than-representative of the
sample of network space our process is attempting to explore.
(In part, this is an indication of why one should not rely
solely on the MPLE). The package contains many algorithmic 
tools to obtain quality inference.
If the sample of networks is
degenerate, our algorithm will fail in its calculation of an
MLE for our data (usually in constructing the Hessian matrix).
%In this case, we can attempt to converge using the Robbins-Monro
%algorithm, minimizing the Euclidean distance between the mean of
%the sample network statistics and the observed network statistics
%as a function of the theta vector.  While this is a fairly
%blunt instrument for estimating the thetas in comparison to
%the MCMC, it still has the potential to get us thetas in the
%right range; then we can run the MCMC again, using our better
%estimates of theta.  A potential problem with this technique is
%that it only works well if the space of sufficient statistics
%is sufficiently regular around our sample network; but as a
%stop-gap measure, it has shown some limited success.
See the references for details, especially Handcock (2003) and Hunter and
Handcock (2006).
}

\references{
Boer, P., Huisman, M., Snijders, T.A.B., and Zeggelink, E.P.H. (2003). 
\emph{StOCNET: an open software system for the advanced statistical analysis 
of social networks.} Version 1.4. Groningen: ProGAMMA / ICS

Handcock, M.S. (2000) \emph{Progress in Statistical Modeling of Drug User and 
  Sexual Networks}, Center for Statistics and the Social Sciences,
     University of Washington.

Handcock, M. S. (2002)
    \emph{Degeneracy and inference for social network models}
    Paper presented at the Sunbelt XXII International Social
Network Conference in New Orleans, LA.

Handcock, M. S. (2003)
    \emph{Assessing Degeneracy in Statistical Models of Social Networks},
    Working Paper \#39, 
Center for Statistics and the Social Sciences,
University of Washington.
\url{www.csss.washington.edu/Papers/wp39.pdf}

Hunter, D. R. and Handcock, M. S. (2006)
    \emph{Inference in curved exponential family models for networks},
   Journal of Computational and Graphical Statistics.
}

\seealso{network, \%v\%, \%n\%, ergm-terms,
         summary.ergm, print.ergm}
\examples{
#
# load the Florentine marriage data matrix
#
data(flo)
#
# attach the sociomatrix for the Florentine marriage data
# This is not yet a network object.
#
flo
#
# Create a network object out of the adjacency matrix
#
flomarriage <- network(flo,directed=FALSE)
flomarriage
#
# print out the sociomatrix for the Florentine marriage data
#
flomarriage[,]
#
# create a vector indicating the wealth of each family (in thousands of lira) 
# and add it as a covariate to the network object
#
flomarriage \%v\% "wealth" <- c(10,36,27,146,55,44,20,8,42,103,48,49,10,48,32,3)
flomarriage
#
# create a plot of the social network
#
plot(flomarriage)
#
# now make the vertex size proportional to their wealth
#
plot(flomarriage, vertex.cex="wealth", main="Marriage Ties")
#
# Use 'data(package = "ergm")' to list the data sets in a
#
data(package="ergm")
#
# Load a network object of the Florentine data
#
data(florentine)
#
# Fit a model where the propensity to form ties between
# families depends on the absolute difference in wealth
#
gest <- ergm(flomarriage ~ edges + absdiff("wealth"))
summary(gest)
#
# add terms for the propensity to form 2-stars and triangles
# of families 
#
gest <- ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle)
summary(gest)

# import synthetic network that looks like a molecule
data(molecule)
# Add a attribute to it to mimic the atomic type
molecule \%v\% "atomic type" <- c(1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,3)
#
# create a plot of the social network
# colored by atomic type
#
plot(molecule, vertex.col="atomic type",vertex.cex=3)

# measure tendency to match within each atomic type
gest <- ergm(molecule ~ edges + kstar(2) + triangle + nodematch("atomic type"),
  MCMCsamplesize=10000)
summary(gest)

# compare it to differential homophily by atomic type
gest <- ergm(molecule ~ edges + kstar(2) + triangle
                              + nodematch("atomic type",diff=TRUE),
  MCMCsamplesize=10000)
summary(gest)
}
\keyword{models}
