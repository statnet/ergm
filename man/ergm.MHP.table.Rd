%  File man/ergm.allstats.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2015 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm.MHP.table}
\alias{ergm.MHP.table}
\title{Table mapping reference,constraints, etc. to Metropolis Hastings Proposals (MHP)}
\description{
		This is a low-level function not intended to be called directly by end users. For information on Metropolis-Hastings proposal methods, \link{ergm_MH_proposals}. This function sets up the table mapping constraints, references, etc. to
		MHproposals.  (It is defined in the scope of local environment)  }
\usage{
 ergm.MHP.table(Class, Reference, Constraints, Priority, Weights, MHP)
	}
\arguments{
  \item{Class}{
   default to "c"
  }
  
  \item{Reference}{
  The reference measure used in the model. For the list of reference measures, see \code{\link{ergm-references}}
  }
	
	\item{Constraints}{
	The constraints used in the model. For the list of constraints, see \code{\link{ergm-constraints}}
  }
	
	\item{Priority}{
	On existence of multiple qualifying MHPs, specifies 
	the priority (-1,0,1) of MHPs to be used.
  }
  
  \item{Weights}{
  The sampling weights on selecting toggles (random, TNT, etc). 
  }
	
	\item{MHP}{
	The matching MHP from the previous arguments. 
	}
}


