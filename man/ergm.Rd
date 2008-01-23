\name{ergm}
\alias{ergm}
\alias{ergm.object}
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
     constraints=~.,
     control=control.ergm(),
     verbose=FALSE, \dots)
}
\arguments{
  \item{formula}{formula; an \R \code{\link{formula}} object, of the form
    \code{y ~ <model terms>},
    where \code{y} is a \code{\link[network]{network}} object or a matrix that can be
    coerced to a \code{\link[network]{network}}  object.  For the details on the possible
    \code{<model terms>}, see \code{\link{ergm-terms}} and Morris, Handcock and
Hunter (2008).  To create a
    \code{\link[network]{network}} object in \R, use the \code{network()} function,
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
  \item{MCMCsamplesize}{count; the number of network statistics,
    randomly drawn from a given distribution on the set of all networks,
    returned by the Metropolis-Hastings algorithm.}
  \item{interval}{count; the number of proposals between sampled statistics.
  }
  \item{maxit}{count; the number of times the parameter for the MCMC
    should be updated by maximizing the MCMC likelihood. At each
    step the parameter is changed to the values that maximizes
    the MCMC likelihood based on the current sample.}
  \item{constraints}{A one-sided formula specifying one or more constraints
    on the support of the distribution of the networks being modeled,
    using syntax similar to the \code{formula} argument. Multiple constraints
    may be given, separated by \dQuote{+} operators.
    Together with the model terms in the formula, the constraints define the distribution of networks being modeled.
        
    It is also possible to specify a proposal function directly by
    passing a string with the function's name. In that case, arguments
    to the proposal should be specified through the \code{prop.args}
    argument to \code{\link{control.ergm}}.
    
    The default is \code{~.}, for an unconstrained model.
    
    The constraint terms currently implemented are
    \describe{
      \item{\code{.} or \code{NULL}}{
	A placeholder for no constraints: all networks of a particular size and type have
	non-zero probability. Cannot be combined with other constraints.
      }
      \item{\code{bd(attribs,maxout,maxin,minout,minin)}}{
	Constrain maximum and minimum vertex degree. See \dQuote{Placing Bounds
    on Degrees} section for more information.
      }
      \item{\code{degrees} and \code{nodedegrees}}{
	Preserve the degree of each vertex of the given network: only networks
	whose vertex degrees are the same as those in the network passed
	in the model formula have non-zero probability.
      }
      \item{\code{degreedist}}{
	Preserve the degree distribution of the given network: only networks
	whose degree distributions are the same as those in the network passed
	in the model formula have non-zero probability.
      }
      \item{\code{indegreedist} and \code{outdegreedist}}{
	Preserve the (respectively) indegree or outdegree distribution of 
        the given network.
      }
      \item{\code{edges}}{
	Preserve the edge count of the given network: only networks
	having the same number of edges as the network passed
	in the model formula have non-zero probability.
      }
%      \item{\code{hamming}}{
%	No clue.
%      }
    }
    Not all combinations of the above are supported.
  }
  \item{control}{A list of control parameters for algorithm
    tuning. Constructed using \code{\link{control.ergm}}. }
  \item{seed}{integer; random number integer seed. Defaults to \code{NULL} to
    use whatever the state of the random number generator is at the time
    of the call.
  }
  \item{verbose}{logical; if this is
    \code{TRUE}, the program will print out additional information,
    including goodness of fit statistics.
  }
  \item{\dots}{Additional
    arguments, to be passed to lower-level functions in the future.
  }
}

\section{Model Terms}{
  The \code{\link{ergm}} function allows the user to explore a large number
  of potential models for their network data.  The
  terms currently supported by the program,
  and a brief description of each is given in the documentation
  \code{\link{ergm-terms}}. 
  In the \code{\link{formula}} for the model, the model terms are various function-like
  calls, some of which require arguments, separated by \code{+} signs.
  For a more detailed understanding of the model terms, see and Morris, Handcock and Hunter (2008).
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
 \item{loglikelihood}{The approximate change in log-likelihood 
 in the last iteration.
  The value is only approximate because it is estimated based 
  on the MCMC random sample.}
 \item{gradient}{The value of the gradient vector of the approximated
  loglikelihood function, evaluated at the maximizer.  This vector
  should be very close to zero.}
 \item{covar}{Approximate covariance matrix for the MLE, based on the inverse
  Hessian of the approximated loglikelihood evaluated at the maximizer.}
 \item{samplesize}{The size of the MCMC sample}
 \item{failure}{Logical:  Did the MCMC estimation fail?}
 \item{mc.se}{MCMC standard error estimates}
 \item{newnetwork}{The final network at the end of the MCMC
 simulation}
 \item{burnin}{If included, the burnin used for the MCMC simulation}
 \item{interval}{If included, the interval used for the MCMC simulation}
 \item{network}{Original network}
 \item{theta.original}{The first value of theta0}
 \item{mplefit}{The MPLE fit as a \code{glm} object.}
 \item{null.deviance}{Deviance of the null model.}
 \item{mle.lik}{The approximate log-likelihood for the MLE.
  The value is only approximate because it is estimated based 
  on the MCMC random sample.}
 \item{etamap}{The set of functions mapping the true parameter theta
 to the canonical parameter eta (irrelevant except in a curved exponential
 family model)}
 \item{degeneracy.value}{Score calculated to assess the degree of 
 degeneracy in the model.}
 \item{degeneracy.type}{Supporting output for \code{degeneracy.value}. Mainly
for internal use.}
 \item{formula}{The original \code{\link{formula}} entered into the \code{\link{ergm}} function.}
 \item{constraints}{Constraints used by original \code{ergm} call}
 \item{prop.weights}{MCMC proposal weights used by original \code{ergm} call
 (part of the \code{\link{control.ergm}} function output).}
 \item{offset}{vector of logical telling which model parameters are to be set
 at a fixed value (i.e., not estimated).}
 \item{drop}{list of terms that were dropped due to extreme values of the
 corresponding statistics on the observed network.}
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
\code{control.ergm(style="Robbins-Monro")} option.
%There seems little reason to implement the Newton-Raphson
%algorithm in C since it moves very quickly in \R and
%it is easier to alter the code if it is written in \R.

The default mechanism for proposing new networks for the MCMC sample 
space is the Metropolis-Hastings algorithm, which simply chooses a dyad 
at random and proposes to toggle that edge; each possible dyad is 
equally likely.  The \code{proposaltype} 
option allow many more complex proposals to be specified.
We have developed and implemented a wide range of algorithms. These are 
described in the documentation for \code{proposaltype}.
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

The package is designed so that the user can add additional 
proposal types. 
}
  
\section{Placing Bounds on Degrees:}{
There are many times when one may wish to condition on the
number of inedges or outedges possessed by a node, either as a
consequence of some intrinsic property of that node (e.g., to control for 
activity or popularity processes), to account
for known outliers of some kind, and thus we wish to limit its indegree, an
intrinsic property of the sampling scheme whence came our data (e.g.,
the survey asked everyone to name only three friends total) or as a
function of the attributes of the nodes to which a node has edges
(e.g., we specify that nodes designated \dQuote{male} have a maximum number
of outdegrees to nodes designated \dQuote{female}).  To accomplish this we
use the \code{constraints} term  \code{bd}.

Let's consider the simple cases first. 
Suppose you want to condition on the total number of degrees
regardless of attributes.  That is, if you had a survey that asked
respondents to name three alters and no more, then you might want to
limit your maximal outdegree to three without regard to any of the
alters' attributes. The argument is then:  

   \code{constraints=~bd(maxout=3)}

Similar calls are used to restrict the number of indegrees
(\code{maxin}), the minimum number of outdegrees
(\code{minout}), and the minimum number of indegrees
(\code{minin}).

You can also set ego specific limits. For example:

   \code{constraints=bd(maxout=rep(c(3,4),c(36,35)))}

limits the first 36 to 3 and the other 35 to 4 outdegrees.

Multiple restrictions can be combined. \code{bd} is very flexible.
In general, the \code{bd} term can contain up to five arguments:

  \preformatted{
    bd(attribs=attribs,
       maxout=maxout,
       maxin=maxin,
       minout=minout,
       minin=minin)
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
required for \dQuote{male} and \dQuote{female} if we want to condition on
the number of ties to both \dQuote{male} and \dQuote{female} partners).
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
\dQuote{sex}), and the attribute sex is represented in ads.sex as
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

% First, you would assign the \dQuote{ego} attribute to all your nodes -
% because every node is an ego, every node gets a \code{TRUE} (or 1) in its
% column in the \code{attribs} matrix.  

% \code{attribs <- matrix(ncol=1,nrow=71,data=c(rep(1,71)))}

% Then you would assign a maximum number of alter out degrees allowed
% to each ego in the same column as in the \code{attribs} matrix.  In this
% case, all egos have a maximum outdegree of 3.

% \code{maxout <- matrix(ncol=1,nrow=71,data=c(rep(3,71)))}

% Then, we can create dummy \code{maxin}, \code{minout},
% \code{minin} matrices (in the future these will be empty
% by default) and create the \code{} argument to pass
% into the \code{\link{ergm}} call.

% To create \dQuote{dummy} \code{maxin}, \code{minout}, and
% \code{minin} matrices to correspond to the \code{attribs} and
% \code{maxout} matrices created above, that impose no limits
% on the amount of edges allowed, our code would look like this:

% \code{maxin <- matrix(NA,ncol=2,nrow=71)
%       minout <- matrix(NA,ncol=2,nrow=71)
%	minin <- matrix(NA,ncol=2,nrow=71)}

%  The argument we pass into the \code{\link{ergm}} call, then would look like this:

%  \code{ = list( attribs=attribs,
%        	         maxout=maxout,
%        	         maxin=maxin,
%        	         minout=minout,
%        	         minin=minin)}

%  This argument passed into \code{\link{ergm}} will now restrict your proposed 
%  networks to those which meet your desired degree conditions.
  

%  \code{ = list( attribs=attribs,
%        	         maxout=maxout,
%        	         maxin=matrix(NA,ncol=1,nrow=71),
%        	         minout=matrix(NA,ncol=1,nrow=71),
%        	         minin=matrix(NA,ncol=1,nrow=71))}

\preformatted{
attribs[n][0] = 1 # just the ego values
maxout[n][0] = minout[n][0] = observed outdegree of n in network
maxin[n][0] = minin[n][0] = observed indegree of n in network
  }
}
 
%\section{Dealing with degeneracy}{
%In order to begin the process of estimating network coefficients,
%we need starting values - guesses at the true values of the
%network statistic coefficients.  The default is to begin the MCMC
%estimation process at the deterministic MPLE values. These
%values are often taken as good-enough final answers by many other
%applications. However recent work has indicated that they are sub-optimal and
%can be dramatically bad.
%
%In using the MPLE values, MCMC MLE often runs into
%problems caused by the inherent instability of the natural parameter
%space of the models (Handcock 2000, 2002, 2003).
%If the initial values for the parameter coefficients
%are off by a very small amount in the wrong direction, the
%result is often a sample of networks that are degenerate -
%that is, networks that are entirely full or entirely empty,
%or that are otherwise less-than-representative of the
%sample of network space our process is attempting to explore.
%(In part, this is an indication of why one should not rely
%solely on the MPLE). The package contains many algorithmic 
%tools to obtain quality inference.
%If the sample of networks is
%degenerate, our algorithm will fail in its calculation of an
%MLE for our data (usually in constructing the Hessian matrix).
%%In this case, we can attempt to converge using the Robbins-Monro
%%algorithm, minimizing the Euclidean distance between the mean of
%%the sample network statistics and the observed network statistics
%%as a function of the theta vector.  While this is a fairly
%%blunt instrument for estimating the thetas in comparison to
%%the MCMC, it still has the potential to get us thetas in the
%%right range; then we can run the MCMC again, using our better
%%estimates of theta.  A potential problem with this technique is
%%that it only works well if the space of sufficient statistics
%%is sufficiently regular around our sample network; but as a
%%stop-gap measure, it has shown some limited success.
%See the references for details, especially Handcock (2003) and Hunter and
%Handcock (2006).
%}

\references{
Boer, P., Huisman, M., Snijders, T.A.B., and Zeggelink, E.P.H. (2003). 
\emph{StOCNET: an open software system for the advanced statistical analysis 
of social networks.} Version 1.4. Groningen: ProGAMMA / ICS

%Handcock, M.S. (2000) \emph{Progress in Statistical Modeling of Drug User and 
%  Sexual Networks}, Center for Statistics and the Social Sciences,
%     University of Washington.

%Handcock, M. S. (2002)
%    \emph{Degeneracy and inference for social network models}
%    Paper presented at the Sunbelt XXII International Social
%Network Conference in New Orleans, LA.

Admiraal R, Handcock MS (2007).
 {\pkg{networksis}: Simulate bipartite graphs with fixed
  marginals through sequential importance sampling}.
 Statnet Project, Seattle, WA.
 Version 1.\newline\url{http://statnetproject.org}.

Bender-deMoll S, Morris M, Moody J (2008).
 {Prototype Packages for Managing and Animating Longitudinal
  Network Data: \pkg{dynamicnetwork} and \pkg{rSoNIA}.}
 {Journal of Statistical Software}, \textbf{24}(7).
 \url{http://www.jstatsoft.org/v24/i07/}.

Boer P, Huisman M, Snijders T, Zeggelink E (2003).
 {StOCNET: an open software system for the advanced statistical
  analysis of social networks}.
 Groningen: ProGAMMA / ICS, version 1.4 edition.

Butts CT (2006).
 {\pkg{netperm}: Permutation Models for Relational Data}.
 Version 0.2, \url{http://erzuli.ss.uci.edu/R.stuff}.

Butts CT (2007).
 {\pkg{sna}: Tools for Social Network Analysis}.
 Version 1.5, \url{http://erzuli.ss.uci.edu/R.stuff}.

Butts CT (2008).
 {\pkg{network}: {A} Package for Managing Relational Data in \R.}
 {Journal of Statistical Software}, \textbf{24}(2).
 \url{http://www.jstatsoft.org/v24/i02/}.

Butts CT, with help~from David~Hunter, Handcock MS (2007).
 {\pkg{network}: Classes for Relational Data}.
 Version 1.2, \url{http://erzuli.ss.uci.edu/R.stuff}.

Goodreau SM, Handcock MS, Hunter DR, Butts CT, Morris M (2008a).
 {A \pkg{statnet} Tutorial.}
 {Journal of Statistical Software}, \textbf{24}(8).
 \url{http://www.jstatsoft.org/v24/i08/}.

Goodreau SM, Kitts J, Morris M (2008{{b}}).
 {Birds of a Feather, or Friend of a Friend? Using Exponential
  Random Graph Models to Investigate Adolescent Social Networks.}
 {Demography}, \textbf{45}, in press.

Handcock, M. S. (2003)
    \emph{Assessing Degeneracy in Statistical Models of Social Networks},
    Working Paper \#39, 
Center for Statistics and the Social Sciences,
University of Washington.
\url{www.csss.washington.edu/Papers/wp39.pdf}

Handcock MS (2003{{b}}).
 {\pkg{degreenet}: Models for Skewed Count Distributions Relevant
  to Networks}.
 Statnet Project, Seattle, WA.
 Version 1.0, \url{http://statnetproject.org}.

Handcock MS, Hunter DR, Butts CT, Goodreau SM, Morris M (2003{{a}}).
 {\pkg{ergm}: {A} Package to Fit, Simulate and Diagnose
  Exponential-Family Models for Networks}.
 Statnet Project, Seattle, WA.
 Version 2, \url{http://statnetproject.org}.

Handcock MS, Hunter DR, Butts CT, Goodreau SM, Morris M (2003{{b}}).
 {\pkg{statnet}: Software Tools for the Statistical Modeling of
  Network Data}.
 Statnet Project, Seattle, WA.
 Version 2, \url{http://statnetproject.org}.

Hunter, D. R. and Handcock, M. S. (2006)
    \emph{Inference in curved exponential family models for networks},
   Journal of Computational and Graphical Statistics.

Hunter DR, Handcock MS, Butts CT, Goodreau SM, Morris M (2008{{b}}).
 {\pkg{ergm}: {A} Package to Fit, Simulate and Diagnose
  Exponential-Family Models for Networks.}
 {Journal of Statistical Software}, \textbf{24}(3).
 \url{http://www.jstatsoft.org/v24/i03/}.

Krivitsky PN, Handcock MS (2007).
 {\pkg{latentnet}: Latent position and cluster models for
  statistical networks}.
 Seattle, WA.
 Version 2, \url{http://statnetproject.org}.

Morris M, Handcock MS, Hunter DR (2008).
 {Specification of Exponential-Family Random Graph Models:
  Terms and Computational Aspects.}
 {Journal of Statistical Software}, \textbf{24}(4).
 \url{http://www.jstatsoft.org/v24/i04/}.
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
