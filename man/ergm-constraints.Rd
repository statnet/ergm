%  File man/ergm-constraints.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2014 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm-constraints}
\alias{ergm-constraints}
\alias{constraints-ergm}
\alias{ergm.constraints}
\alias{constraints.ergm}

\alias{ConstraintImplications}
\alias{InitConstraint.bd}
\alias{InitConstraint.degreedist}
\alias{InitConstraint.degrees}
\alias{InitConstraint.odegrees}
\alias{InitConstraint.idegrees}
\alias{InitConstraint.b1degrees}
\alias{InitConstraint.b2degrees}
\alias{InitConstraint.blockdiag}
\alias{InitConstraint.edges}
\alias{InitConstraint.hamming}
\alias{InitConstraint.idegreedist}
\alias{InitConstraint.nodedegrees}
\alias{InitConstraint.observed}
\alias{InitConstraint.odegreedist}
\alias{InitConstraint.fixedas}
\alias{InitConstraint.fixallbut}

\docType{methods}

\title{Sample Space Constraints for Exponential Family Random Graph Models}
\description{
  \code{\link{ergm}} is used to fit linear exponential random graph
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
\section{Constraints implemented in  the \code{\link[=ergm-package]{ergm}} package}{
  \describe{
    \item{\code{.} or \code{NULL}}{
      A placeholder for no constraints: all networks of a
      particular size and type have non-zero probability.
      Cannot be combined with other constraints.
    }
    \item{\code{bd(attribs,maxout,maxin,minout,minin)}}{
      Constrain maximum and minimum vertex degree. See
      \dQuote{Placing Bounds
        on Degrees} section for more information.
    }
    \item{\code{blockdiag(attrname)}}{
      Force a block-diagonal structure on the network. Only dyads
      \eqn{(i,j)} for which \code{attrname(i)==attrname(j)} can have edges.
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
    \item{\code{idegreedist} and \code{odegreedist}}{
      Preserve the (respectively) indegree or outdegree distribution of 
      the given network.
    }
    \item{\code{edges}}{
      Preserve the edge count of the given network: only networks
      having the same number of edges as the network passed
      in the model formula have non-zero probability.
    }
    \item{\code{observed}}{
      Preserve the observed dyads of the given network.
    }
    \item{\code{fixedas(present,absent)}}{
      Preserve the edges in 'present' and preclude the edges in 'absent'. 
      Both 'present' and 'absent' can take input object as edgelist and 
      network, the latter will convert to the corresponding edgelist.
 	}
 	
 	\item{\code{fixallbut(free.dyads)}}{
      Preserve the dyad status in all but free.dyads. free.dyads can 
      take input object as edgelist and network, the latter will 
      convert to the corresponding edgelist.
 	}
  Not all combinations of the above are supported.
  }
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
 \url{http://www.jstatsoft.org/v24/i08/}.

Hunter, D. R. and Handcock, M. S. (2006)
    \emph{Inference in curved exponential family models for networks},
   Journal of Computational and Graphical Statistics.

Hunter DR, Handcock MS, Butts CT, Goodreau SM, Morris M (2008b).
 \pkg{ergm}: A Package to Fit, Simulate and Diagnose
  Exponential-Family Models for Networks.
 \emph{Journal of Statistical Software}, 24(3).
 \url{http://www.jstatsoft.org/v24/i03/}.

Krivitsky PN (2012). Exponential-Family Random Graph Models for Valued
Networks. \emph{Electronic Journal of Statistics}, 2012, 6,
1100-1128. \href{http://dx.doi.org/10.1214/12-EJS696}{\code{doi:10.1214/12-EJS696}}
 
Morris M, Handcock MS, Hunter DR (2008).
 Specification of Exponential-Family Random Graph Models:
  Terms and Computational Aspects.
 \emph{Journal of Statistical Software}, 24(4).
 \url{http://www.jstatsoft.org/v24/i04/}.
  
}

\keyword{models}
