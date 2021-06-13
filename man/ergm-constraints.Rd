%  File man/ergm-constraints.Rd in package ergm, part of the
%  Statnet suite of packages for network analysis, https://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  https://statnet.org/attribution .
%
%  Copyright 2003-2021 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm-constraints}
\alias{ergm-constraints}
\alias{constraints-ergm}
\alias{ergm.constraints}
\alias{constraints.ergm}

\alias{Dyads}
\alias{bd}
\alias{degreedist-constraint}
\alias{degrees}
\alias{dyadnoise}
\alias{odegrees}
\alias{idegrees}
\alias{b1degrees}
\alias{b2degrees}
\alias{edges-constraint}
\alias{hamming-constraint}
\alias{idegreedist}
\alias{nodedegrees}
\alias{observed}
\alias{odegreedist}
\alias{fixedas}
\alias{fixallbut}
\alias{egocentric}
\alias{blocks}

\docType{methods}

\title{Sample Space Constraints for Exponential-Family Random Graph Models}
\description{
  \code{\link{ergm}} is used to fit exponential-family random graph
  models (ERGMs), in which
  the probability of a given network, \eqn{y}, on a set of nodes is 
  \eqn{h(y) \exp\{\eta(\theta) \cdot
    g(y)\}/c(\theta)}, where
  \eqn{h(y)} is the reference measure (usually \eqn{h(y)=1}),
  \eqn{g(y)} is a vector of network statistics for \eqn{y},
  \eqn{\eta(\theta)} is a natural parameter vector of the same 
  length (with \eqn{\eta(\theta)=\theta} for most terms), and \eqn{c(\theta)} is the
  normalizing constant for the distribution.

  This page describes the constraints (the networks \eqn{y} for which \eqn{h(y)>0}) 
  that are included with the
  \code{\link[=ergm-package]{ergm}} package. Other packages may add new
  constraints.

}
\section{Constraints formula}{
  
  A constraints formula is a one- or two-sided formula whose left-hand
  side is an optional direct selection of the \code{InitErgmProposal} function
  and whose right-hand
  side is a series of one or more terms separated by \dQuote{\code{+}}
  and \dQuote{\code{-}} operators, specifying the constraint.

  The sample space (over and above the reference distribution) is
  determined by iterating over the constraints terms from left to right,
  each term updating it as follows:
  \itemize{
    \item If the constraint introduces complex dependence structure
    (e.g., constrains degree or number of edges in the network), then
    this constraint always restricts the sample space. It may only have
    a \dQuote{\code{+}} sign.
    
    \item If the constraint only restricts the set of dyads that may
    vary in the sample space (e.g., block-diagonal structure or fixing
    specific dyads at specific values) and has a \dQuote{\code{+}} sign,
    the set of dyads that may vary is restricted to those that may vary
    according to this constraint \emph{and} all the constraints to date.

    \item If the constraint only restricts the set of dyads that may
    vary in the sample space but has a \dQuote{\code{-}} sign,
    the set of dyads that may vary is expanded to those that may vary
    according to this constraint \emph{or} all the constraints up to date.    
  }

  For example, a constraints formula \code{~a-b+c-d} with all
  constraints dyadic will allow dyads permitted by either `a` or `b` but
  only if they are also permitted by `c`; as well as all dyads permitted
  by `d`. If `A`, `B`, `C`, and `D` were logical matrices, the matrix of
  variable dyads would be equal to `((A|B)&C)|D`.
  
  Terms with a positive sign can be viewed as \dQuote{adding} a
  constraint while those with a negative sign can be viewed as
  \dQuote{relaxing} a constraint.
}

\section{Constraints implemented in  the \code{\link[=ergm-package]{ergm}} package}{
  \describe{
    \item{\code{.} or \code{NULL} (dyad-independent)}{
      A placeholder for no constraints: all networks of a
      particular size and type have non-zero probability.
      Cannot be combined with other constraints.
    }
    \item{\code{Dyads(fix=NULL, vary=NULL)} (dyad-independent)}{
      This is an \dQuote{operator} constraint that takes one or two \code{\link{ergm}} formulas. These formulas should contaion only dyad-independent terms. For the terms in the \code{fix=} formula, dyads that affect the network statistic (i.e., have nonzero change statistic) for \emph{any} the terms will be fixed at their current values. For the terms in the \code{vary=} formula, only those that change \emph{at least one} of the terms will be allowed to vary, and all others will be fixed. If both formulas are given, the dyads that vary either for one or for the other will be allowed to vary. Note that a formula passed to \code{Dyads} without an argument name will default to \code{fix=}.
    }
    \item{\code{bd(attribs,maxout,maxin,minout,minin)}}{
      Constrain maximum and minimum vertex degree. See
      \dQuote{Placing Bounds
        on Degrees} section for more information.
    }
    \item{\code{degrees} and \code{nodedegrees}}{
      Preserve the degree of each vertex of the given network: only networks
      whose vertex degrees are the same as those in the network passed
      in the model formula have non-zero probability. If the network is
      directed, both indegree and outdegree are preserved.
    }
    \item{\code{odegrees}, \code{idegrees}, \code{b1degrees}, \code{b2degrees}}{
      For directed networks, \code{odegrees} preserves the outdegree of each vertex of the given
      network, while allowing indegree to vary, and conversely for
      \code{idegrees}. \code{b1degrees} and \code{b2degrees} perform a
      similar function for bipartite networks.
    }
    \item{\code{degreedist}}{
      Preserve the degree distribution of the given network: only networks
      whose degree distributions are the same as those in the network passed
      in the model formula have non-zero probability.
    }
    \item{\code{dyadnoise(p01,p10)}}{
      A soft constraint to adjust the sampled distribution for
      dyad-level noise with known perturbation probabilities. It is
      assumed that the observed LHS network is a noisy observation of
      some unobserved true network, with \code{p01} giving the dyadwise
      probability of erroneously observing a tie where the true network
      had a non-tie and \code{p10} giving the dyadwise probability of
      erroneously observing a nontie where the true network had a tie.

      \code{p01} and \code{p10} can be either both be scalars or or both
      be adjacency matrices of the same dimension as that of the LHS
      network giving these probabilities.

      See Karwa et al. (2016) for an application.
    }
    \item{\code{idegreedist} and \code{odegreedist}}{
      Preserve the (respectively) indegree or outdegree distribution of 
      the given network.
    }
    \item{\code{edges}}{
      Preserve the edge count of the given network: only networks
      having the same number of edges as the network passed
      in the model formula have non-zero probability.
    }
    \item{\code{observed} (dyad-independent)}{
      Preserve the observed dyads of the given network.
    }
    \item{\code{fixedas(present,absent)} (dyad-independent)}{
      Preserve the edges in 'present' and preclude the edges in 'absent'. 
      Both 'present' and 'absent' can take input object as edgelist and 
      network, the latter will convert to the corresponding edgelist.
 	}
 	
 	\item{\code{fixallbut(free.dyads)} (dyad-independent)}{
      Preserve the dyad status in all but free.dyads. free.dyads can 
      take input object as edgelist and network, the latter will 
      convert to the corresponding edgelist.
    }

    \item{\code{egocentric(attr = NULL, direction = c("both","out","in"))} (dyad-independent)}{
      Preserve values of dyads incident on vertices with attribute \code{attr} (see \link[=nodal_attributes]{Specifying Vertex Attributes and Levels} for details) being \code{TRUE} or if \code{attrname} is \code{NULL}, the vertex attribute \code{"na"} being `FALSE`. For directed networks, \code{direction=="out"} only preserves the out-dyads of those actors, and \code{direction=="in"} preserves their in-dyads.
    }

    \item{\code{blocks(attr = NULL, levels = NULL, levels2 = FALSE, b1levels = NULL, b2levels = NULL)} (dyad-independent)}{
      Constrain "blocks" of dyads; any dyad whose toggle would produce a nonzero change statistic for a \code{nodemix} term with the same arguments will be fixed.  Note that the \code{levels2} argument has a different default value for \code{blocks} than it does for \code{nodemix}.
    }

    \item{\code{blockdiag(attr)} (dyad-independent)}{
      Force a block-diagonal structure (and its bipartite analogue) on
      the network. Only dyads \eqn{(i,j)} for which
      \code{attr(i)==attr(j)} can have edges. See \link[=nodal_attributes]{Specifying Vertex attributes and Levels} (\verb{? nodal_attributes}) for the ways to specify nodal attributes and expressions.

      Note that the current implementation requires that blocks be
      contiguous for \dQuote{unipartite} graphs, and for bipartite
      graphs, they must be contiguous within a partition and must have
      the same ordering in both partitions. (They do not, however,
      require that all blocks be represented in both partitions, but
      those that overlap must have the same order.)

      If multiple block-diagonal constraints are given, or if
      \code{attr} is a vector with multiple attribute names, blocks
      will be constructed on \emph{all} attributes matching.
    }
  }
  Not all combinations of the above are supported.
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

%\preformatted{
%attribs[n][0] = 1 # just the ego values
%maxout[n][0] = minout[n][0] = observed outdegree of n in network
%maxin[n][0] = minin[n][0] = observed indegree of n in network
%  }
}
 

\references{

Goodreau SM, Handcock MS, Hunter DR, Butts CT, Morris M (2008a).
 A \pkg{statnet} Tutorial.
 \emph{Journal of Statistical Software}, 24(8).
 \url{https://www.jstatsoft.org/v24/i08/}.

Hunter, D. R. and Handcock, M. S. (2006)
    \emph{Inference in curved exponential family models for networks},
   Journal of Computational and Graphical Statistics.

Hunter DR, Handcock MS, Butts CT, Goodreau SM, Morris M (2008b).
 \pkg{ergm}: A Package to Fit, Simulate and Diagnose
  Exponential-Family Models for Networks.
 \emph{Journal of Statistical Software}, 24(3).
 \url{https://www.jstatsoft.org/v24/i03/}.

Karwa V, Krivitsky PN, and Slavkovi\'{c} AB (2016). Sharing Social Network
 Data: Differentially Private Estimation of Exponential-Family Random
 Graph Models. \emph{Journal of the Royal Statistical Society, Series
 C}, 66(3): 481-500. \doi{10.1111/rssc.12185}
 
Krivitsky PN (2012). Exponential-Family Random Graph Models for Valued
Networks. \emph{Electronic Journal of Statistics}, 6,
1100-1128. \doi{10.1214/12-EJS696}
 
Morris M, Handcock MS, Hunter DR (2008).
 Specification of Exponential-Family Random Graph Models:
  Terms and Computational Aspects.
 \emph{Journal of Statistical Software}, 24(4).
 \url{https://www.jstatsoft.org/v24/i04/}.
  
}

\keyword{models}
