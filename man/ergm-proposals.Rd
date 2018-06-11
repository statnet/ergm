%  File man/ergm-proposals.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2018 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm-proposals}
\alias{ergm-proposals}
\alias{InitErgmProposal}
\alias{InitErgmProposal.CondDegree}
\alias{InitErgmProposal.CondDegreeDist}
\alias{InitErgmProposal.CondInDegreeDist}
\alias{InitErgmProposal.CondOutDegreeDist}
\alias{InitErgmProposal.ConstantEdges}
\alias{InitErgmProposal.TNT}
\alias{InitErgmProposal.randomtoggle}
\alias{InitErgmProposal.randomtoggleNonObserved}
\alias{InitErgmProposal.NonObservedTNT}
\alias{InitErgmProposal.CondB1Degree}
\alias{InitErgmProposal.CondB2Degree}
\alias{InitErgmProposal.CondInDegree}
\alias{InitErgmProposal.CondOutDegree}
\alias{InitErgmProposal.blockdiag}
\alias{InitErgmProposal.blockdiagNonObserved}
\alias{InitErgmProposal.blockdiagTNT}
\alias{InitErgmProposal.CondDegreeMix}
\alias{InitErgmProposal.fixedas}
\alias{InitErgmProposal.fixallbut}
\alias{InitErgmProposal.HammingConstantEdges}
\alias{InitErgmProposal.HammingTNT}
\alias{InitErgmProposal.nobetweengroupties}

\alias{NonObservedTNT}
\alias{InitErgmProposal.blockdiagNonObservedTNT}
\alias{InitErgmProposal.fixallbutTNT} 
\alias{InitErgmProposal.fixedasTNT}


\docType{methods}

\title{Metropolis-Hastings Proposal Methods for ERGM MCMC}
\description{
  \code{\link{ergm}}  uses a Metropolis-Hastings (MH) algorithm to control the behavior of the Markov Chain Monte Carlo (MCMC) for sampling networks.  The MCMC chain is intended to step around the sample space of possible networks, selecting a network at regular intervals to evaluate the statistics in the model.  For each MCMC step, \eqn{n} (\eqn{n=1} in the simple case) toggles are proposed to change the dyad(s) to the opposite value. The probability of accepting the proposed change is determined by the MH acceptance ratio.  The role of the different MH methods implemented in \code{\link{ergm}} is to vary how the sets of dyads are selected for toggle proposals.  This is used in some cases to improve the performance (speed and mixing) of the algorithm, and in other cases to constrain the sample space.
}

\section{MH proposal methods implemented in the \code{\link[=ergm-package]{ergm}} package}{

\bold{MH proposals for non-constrained ergm models}

  \describe{
  
        \item{InitErgmProposal.randomtoggle}{Propose a randomly selected dyad to toggle.
        }
        
        \item{InitErgmProposal.TNT}{Default MH algorithm. Stratifies the population of dyads
        	by edge status: those having ties and those having no ties (hence T/NT).  
        	This is useful for improving performance in sparse networks,
        	because it gives at least 50\% chance of proposing a toggle of an existing edge.
        }
             
}

\bold{MH proposals for constrained ergm models}


    \describe{

         \item{InitErgmProposal.blockdiag}{MHp for \eqn{constraints= ~blockdiag}. Select a diagonal block according to the weight, then randomly select a dayd within the block for the toggle proposal.
        }
        \item{InitErgmProposal.blockdiagNonObserved}{MHp for \eqn{constraints= ~blockdiag+observed}. Similar to InitErgmProposal.blockdiag, but
        applied only to missing dyads.
        }
        
         \item{InitErgmProposal.blockdiagNonObservedTNT}{Similar to InitErgmProposal.blockdiagNonObserved, except that it selects ties and non-ties for proposed toggles (in the block by construction) with equal probability.  Like the unconstrained TNT proposal, this is useful for improving performance in sparse networks.
        }
        
        
        \item{InitErgmProposal.blockdiagTNT}{MHp for \eqn{constraints= ~blockdiag}. Similar to InitErgmProposal.blockdiag, except that it selects ties
        and non-ties for proposed toggles (in the block by construction) with equal probability.  Like the unconstrained TNT proposal, this is useful for improving performance in sparse networks.
        }
        
         \item{InitErgmProposal.CondB1Degree}{MHp for \eqn{constraints= ~b1degrees}. For bipartite networks, randomly select an edge {B1i, B2j} and an empty dyad with the same node B1i, {B1i, B2k}, and propose to toggle both {B1i, B2j} and {B1i, B2k}. This ensures that the degrees of individual nodes in mode 1 are preserved.
        }
        \item{InitErgmProposal.CondB2Degree}{MHp for \eqn{constraints= ~b2degrees}. For bipartite network, randomly select an edge {B1j, B2i} and an empty dyad with the same node B2i, {B1k, B2i}, and propose to toggle both {B1j, B2i} and {B1k, B2i}. This ensures that the degrees of individual nodes in mode 2 are preserved.
        }
        
        \item{InitErgmProposal.CondDegree}{MHp for \eqn{constraints= ~degree}. Propose either 4 toggles (MH_CondDegreeTetrad) or 6 toggles (MH_CondDegreeHexad) at once. For undirected networks, propose 4 toggles (MH_CondDegreeTetrad). MH_CondDegreeTetrad selects two edges with no nodes in common, A1-A2 and B1-B2, s.t. A1-B2 and B1-A2 are not edges, and propose to replace the former two by the latter two. MH_CondDegreeHexad selects three edges A1->A2, B1->B2, C1->C2 at random and rotate them to A1->B2, B1->C2, and C1->A2.
        }
        \item{InitErgmProposal.CondDegreeDist}{MHp for \eqn{constraints= ~degreedist}. Randomly select a node (T) and its edge (E).  If the head node of the edge (H) has 1 degree more than another randomly select node (A), and A is disconnected to both T and H, then propose to toggle E and the dyad between T and A.
        }
        
        \item{InitErgmProposal.CondDegreeMix}{MHp for \eqn{constraints= ~degreesmix}. Similar to \code{InitErgmProposal.CondDegree}, except that the toggle is proposed only if the mixing matrix of degrees is preserved before and after the toggle.}
        
        \item{InitErgmProposal.ConstantEdges}{MHp for \eqn{constraints= ~edges}. Propose pairs of toggles that keep number of edges the same. This is done by (a) choosing an existing edge at random; (b) repeatedly choosing dyads at random until one is found that does not have an edge; and (c) proposing toggling both these dyads. Note that step (b) will be very inefficient if the network is nearly complete, so this proposal is NOT recommended for such networks. However, most network datasets are sparse, so this is not likely to be an issue.
        }   
        
        \item{InitErgmProposal.CondInDegreeDist}{MHp for \eqn{constraints= ~idegreedist}. For directed networks, similar to InitErgmProposal.CondDegreeDist, except for indegree case
        }
        \item{InitErgmProposal.CondOutDegreeDist}{MHp for \eqn{constraints= ~odegreedist}. For directed networks, similar to InitErgmProposal.CondDegreeDist, except for outdegree case
        }

	\item{InitErgmProposal.fixedas}{MHp for \eqn{constraints=~fixedas(present,absent)}. Select a random dyad that is not in either 'present' edgelist or 'absent' edgelist to toggle. Edges in 'present' and empty dyads in 'absent' are remained fixed.}
	
	\item{InitErgmProposal.fixedasTNT}{Similar to InitErgmProposal.fixedas, except that it selects ties
        and non-ties for proposed toggles with equal probability.  Like the unconstrained TNT proposal, this is useful for improving performance in sparse networks.}	
	
			
	\item{InitErgmProposal.fixallbut}{MHp for \eqn{constraints=~fixallbut(free.dyads)}. Select a random dyad that is in free.dyads edgelist to toggle.}
  
  \item{InitErgmProposal.fixallbutTNT}{Similar to InitErgmProposal.fixallbut, except that it selects ties and non-ties for proposed toggles with equal probability.  Like the unconstrained TNT proposal, this is useful for improving performance in sparse networks.}
  
  
        \item{InitErgmProposal.randomtoggleNonObserved}{MHp for \eqn{constraints= ~observed}. Randomly select a missing/non-observed dyad and propose a toggle.
        }

	\item{InitErgmProposal.NonObservedTNT}{Similar to InitErgmProposal.randomtoggleNonObserved, except that it selects ties and non-ties for proposed toggles with equal probability.  Like the unconstrained TNT proposal, this is useful for improving performance in sparse networks.
        }
	
        \item{InitErgmProposal.CondInDegree}{MHp for \eqn{constraints= ~idegrees}. For directed networks, randomly select two dyads with a common head node, one having an edge one not, and propose to swap the tie from one tail to the other.
        }
        \item{InitErgmProposal.CondOutDegree}{MHp for \eqn{constraints= ~odegrees}. For directed networks, randomly select two dyads with a common tail node, one having an edge and one not, and propose to swap the tie from one head to the other.
        }

       }}

\seealso{\code{\link[=ergm-package]{ergm}} package, \code{\link{ergm}}, \code{\link{ergm-constraints}},  \code{\link{ergm_proposal}}}

  
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
1100-1128. \doi{10.1214/12-EJS696}
 
Morris M, Handcock MS, Hunter DR (2008).
 Specification of Exponential-Family Random Graph Models:
  Terms and Computational Aspects.
 \emph{Journal of Statistical Software}, 24(4).
 \url{http://www.jstatsoft.org/v24/i04/}.
  
}

\keyword{models}



