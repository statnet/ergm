\name{ergm-constraints}

\alias{InitMHP.CondDegree}
\alias{InitMHP.CondDegreeDist}
\alias{InitMHP.CondInDegreeDist}
\alias{InitMHP.CondOutDegreeDist}
\alias{InitMHP.ConstantEdges}
%\alias{InitMHP.HammingConstantEdges}
%\alias{InitMHP.HammingTNT}
\alias{InitMHP.TNT}
\alias{InitMHP.TNT10}
%\alias{InitMHP.nobetweengroupties}
\alias{InitMHP.randomtoggle}
\alias{InitMHP.randomtoggleNonObserved}
\alias{InitMHP.CondB1Degree}
\alias{InitMHP.CondB2Degree}
\alias{InitMHP.CondInDegree}
\alias{InitMHP.CondOutDegree}
\alias{InitMHP.blockdiag}
\alias{InitMHP.blockdiagNonObserved}
\alias{InitMHP.blockdiagTNT}
%\alias{InitWtMHP.StdNormal}
\alias{InitMHP.CondDegreeMix}
\alias{InitMHP.TriNT}

\docType{methods}

\title{Metropolis-Hastings Proposal Methods for ERGM MCMC}
\description{
  \code{\link{ergm}}  uses Metropolis-Hastings (MH) algorithm in Markov Chain Monte Carlo (MCMC) for sampling networks. In each MCMC step, \eqn{n} (\eqn{n=1} in most cases) toggles are proposed to change the dyads to their opposite values. Whether the proposed changes are accepted or not are decided by the MH acceptance ratio. Different MH methods decide which sets of dyads are selected to propose toggles.
}

\section{MH proposal methods implemented in the \code{\link[=ergm-package]{ergm}} package}{

\bold{MH proposals for non-constrained ergm models}

  \describe{
  
        \item{InitMHP.randomtoggle}{Propose a randomly selected dyad to toggle.
        }
        
        \item{InitMHP.TNT}{Default MH algorithm. Gives at least 50\% chance of
   proposing a toggle of an existing edge, as opposed to simple random toggles that rarely do so in sparse networks (comments from code)
        }
             
        \item{InitMHP.TNT10}{(Experimental) Propose 10 TNT steps at once, but this seems flawed currently because it does not correctly update network quantities like nedges after each of the 10 proposed toggles. (comments from code)
        }
        
        \item{InitMHP.TriNT}{Similar to \code{InitMHP.TNT}, except that once the toggle been proposed, randomly select a third node and gives \eqn{50\%} of chance to toggle each dyad that connects the third node to the two nodes from the previous proposed toggle.}
        

}

\bold{MH proposals for constrained ergm models}


    \describe{

         \item{InitMHP.blockdiag}{MHp for \eqn{constraints= ~blockdiag}. Select a diagonal block according to the weight, then randomly select a dayd and propose to toggle.
        }
        \item{InitMHP.blockdiagNonObserved}{MHp for \eqn{constraints= ~blockdiag+observed}. Similar to InitMHP.blockdiag, except that only among missing dyads.
        }
        \item{InitMHP.blockdiagTNT}{MHp for \eqn{constraints= ~blockdiag}. Similar to InitMHP.blockdiag, except that it proposes a toggle of an existing edge (in the block by construction) at least 50\% chance,
        }
        
         \item{InitMHP.CondB1Degree}{MHp for \eqn{constraints= ~b1degrees}. For bipartite network, randomly select an edge A1->A2, a different node in the second mode, B2, and propose to toggle A1->A2 and toggle A1->B2 at once. No matter which mode A1,A2 belong to, the degrees of nodes in mode 1 are preserved  (comments from code)
        }
        \item{InitMHP.CondB2Degree}{MHp for \eqn{constraints= ~b2degrees}. For bipartite network, randomly select an edge B1->B2, a different node in the second mode, A2, and propose to toggle B1->B2 and toggle B1->A2 at once. No matter which mode B1,B2 belong to, the degrees of nodes in mode 2 are preserved
        }
        
        \item{InitMHP.CondDegree}{MHp for \eqn{constraints= ~degree}. Propose either 4 toggles (MH_CondDegreeTetrad) or 6 toggles (MH_CondDegreeHexad) at once. For undirected network, propose 4 toggles (MH_CondDegreeTetrad) MH_CondDegreeTetrad: Select two edges with no nodes in common, A1-A2 and B1-B2, s.t. A1-B2 and B1-A2 are not edges, and propose to replace the former two by the latter two. MH_CondDegreeHexad: Select three edges A1->A2, B1->B2, C1->C2 at random and rotate them to A1->B2, B1->C2, and C1->A2.
        }
        \item{InitMHP.CondDegreeDist}{MHp for \eqn{constraints= ~degreedist}. Randomly select a node (T) and its edge (E), if the head node of the edge (H) has 1 degree more than another randomly select node (A), and A is disconnected to both T and H, then propose to toggle E and the dyad between T and A.
        }
        
        \item{InitMHP.CondDegreeMix}{MHp for \eqn{constraints= ~degreesmix}. Similar to \code{InitMHP.CondDegree}, except that propose the toggle only if the mixing matrix of degrees are preserved before and after the toggle.}
        
        \item{InitMHP.ConstantEdges}{MHp for \eqn{constraints= ~edges}. Propose pairs of toggles that keep number of edges the same. This is done by (a) choosing an existing edge at random; (b) repeatedly choosing dyads at random until one is found that does not have an edge; and (c) proposing toggling both these dyads. Note that step (b) will be very inefficient if the network is nearly complete, so this proposal is NOT recommended for such networks. However, most network datasets are sparse, so this is not likely to be an issue.
        }   
        
        \item{InitMHP.CondInDegreeDist}{MHp for \eqn{constraints= ~idegreedist}. For directed network, similar to InitMHP.CondDegreeDist, except for indegree case
        }
        \item{InitMHP.CondOutDegreeDist}{MHp for \eqn{constraints= ~odegreedist}. For directed network, similar to InitMHP.CondDegreeDist, except for outdegree case
        }

%        \item{InitMHP.HammingConstantEdges}{}
%        \item{InitMHP.HammingTNT}{}

%        \item{InitMHP.nobetweengroupties}{}

        
        \item{InitMHP.randomtoggleNonObserved}{MHp for \eqn{constraints= ~observed}. Randomly select a missing/non-observed dyad and propose a toggle.
        }
       
        \item{InitMHP.CondInDegree}{MHp for \eqn{constraints= ~idegrees}. For directed network, randomly select two dyads with a common head node, one has an edge one doesn't, and propose to swap the tie from one tail to the other.
        }
        \item{InitMHP.CondOutDegree}{MHp for \eqn{constraints= ~odegrees}. For directed network, randomly select two dyads with a common tail node, one has an edge one doesn't, and propose to switch the tie from one head to the other.
        }

%        \item{InitWtMHP.StdNormal}{}

       }}

\seealso{\code{\link[=ergm-package]{ergm}} package, \code{\link{ergm}}, \code{\link{ergm-constraints}}}

  
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



