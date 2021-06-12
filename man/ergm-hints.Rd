%  File man/ergm-hints.Rd in package ergm, part of the
%  Statnet suite of packages for network analysis, https://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  https://statnet.org/attribution .
%
%  Copyright 2003-2021 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm-hints}
\alias{ergm-hints}
\alias{hints}
\alias{sparse}
\alias{strat}

\docType{methods}

\title{Hints for Exponential-Family Random Graph Models}

\description{
  These hints may be used to control proposal probabilities.
}

\section{Hints implemented in the \code{\link[=ergm-package]{ergm}} package}{
  \describe{
    \item{\code{sparse} (dyad-independent)}{
      The network is sparse. This typically results in a Tie-Non-Tie (TNT) proposal regime.
    }
    \item{\code{strat(attr=NULL, pmat=NULL, empirical=FALSE)} (dyad-independent)}{
      The dyads in the network are stratified according to an attribute combination.

      This typically results in stratifying proposals by mixing type on a vertex attribute.

      Specifically, the user may pass a vertex attribute \code{attr} as an argument 
      (the default for \code{attr} gives every vertex the same attribute value),
      and may also pass a matrix of weights \code{pmat} (the default for 
      \code{pmat} gives equal weight to each mixing type).  See 
      \link[=nodal_attributes]{Specifying Vertex Attributes and Levels} for details on 
      specifying vertex attributes.  The matrix \code{pmat}, if specified, must have
      the same dimensions as a mixing matrix for the network and attribute under 
      consideration, and the correspondence between rows and columns of \code{pmat}
      and values of \code{attr} is the same as for a mixing matrix.
                  
      The interpretation is that \code{pmat[i,j]/sum(pmat)} is the probability of 
      proposing a toggle for mixing type \code{(i,j)}.  (For undirected, unipartite
      networks, \code{pmat} is first symmetrized, and then entries below the diagonal 
      are set to zero.  Only entries on or above the diagonal of the symmetrized 
      \code{pmat} are considered when making proposals.  This accounts for the 
      convention that mixing is undirected in an undirected, unipartite network: 
      a tail of type \code{i} and a head of type \code{j} has the same mixing type 
      as a tail of type \code{j} and a head of type \code{i}.)
      
      As an alternative way of specifying \code{pmat}, the user may pass 
      \code{empirical=TRUE} to use the mixing matrix of the network beginning 
      the MCMC chain as \code{pmat}.  In order for this to work,
      that network should have a reasonable (in particular, nonempty) edge set.
      
      While some mixing types may be assigned zero proposal probability
      (either with a direct specification of \code{pmat} or with \code{empirical=TRUE}),
      this will not be recognized as a constraint by all components of \code{ergm},
      and should be used with caution.
    }
  }
}
