%  File man/san.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2017 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{MHproposal}

\alias{MHproposal}
\alias{MHproposal.MHproposal}
\alias{MHproposal.NULL}
\alias{MHproposal.character}
\alias{MHproposal.ergm}
\alias{MHproposal.formula}

\title{Functions to initialize the MHproposal object}
\usage{
                      
\method{MHproposal}{character}(object, 
                   arguments, 
                   nw, \dots,
                   response=NULL, 
                   reference=reference
                   )  
                      
                      
\method{MHproposal}{formula}(object, 
                   arguments, 
                   nw, 
                   weights="default", 
                   class="c", 
                   reference=~Bernoulli, 
                   response=NULL,
                  \dots)

\method{MHproposal}{ergm}(object, \dots,
                   constraints=NULL, 
                   arguments=NULL, 
                   nw=NULL, 
                   weights=NULL,
                   class="c", 
                   reference=NULL, 
                   response=NULL)
}
\description{
S3 Functions that initialize the Metropolis-Hastings Proposal (MHproposal) object using the InitMHP.* function that corresponds to the name given in 'object'.  These functions are not generally called directly by the user.  See \link{ergm_MH_proposals} for general explanation and lists of available Metropolis-Hastings proposal types. 
}
\arguments{
  \item{object}{Either a character, a \code{\link{formula}} or an \code{\link{ergm}} object.
    The \code{\link{formula}} should be of the form \code{y ~ <model terms>},
    where \code{y} is a network object or a matrix that can be
    coerced to a \code{\link[network]{network}} object.  
  }
  \item{nw}{The network object originally given to \code{\link{ergm}} via 'formula'
  }
  
  \item{weights}{Specifies the method used to allocate probabilities of being proposed to dyads; options are "TNT", "TNT10", "random", "nonobserved" and "default"; default="default"}
  
  \item{arguments}{A list of parameters used by the Init.MHP routines}
  
  \item{response}{\emph{EXPERIMENTAL.} Name of the edge attribute whose value is to be
    modeled. Defaults to \code{NULL} for simple presence or absence. }

  \item{reference}{\emph{EXPERIMENTAL.} One-sided formula whose RHS
    gives the reference measure to
    be used. (Defaults to \code{~Bernoulli}.)}

	\item{class}{The class of the proposal; choices include "c", "f", and "d" default="c".
	}
  \item{constraints}{A one-sided formula specifying one or more constraints
    on the support of the distribution of the networks being
    simulated. See the documentation for a similar argument for
    \code{\link{ergm}}  and see \link[=ergm-constraints]{list of
    implemented constraints} for more information.
  }
  \item{\dots}{Further arguments passed to other functions.
  }
}

%\details{}

\value{
   Returns an MHproposal object: a list with class \code{'MHProposal'} containing the following named elements:
   \itemize{
     \item name   : the C name of the proposal
     \item{inputs : NULL (I think - the only non-null value returned by the InitMH
              is for <nobetweengroupties>, but this isn't included in the 
              look-up table}
     \item package: shared library name where the proposal can be found (usually "ergm")
     \item{arguments: list of arguments passed to the InitMHP function; in particular,
       \itemize{
        \item constraints: list of constraints
        \item{constraints$bd: the list of parameters to bound degree in the fitting process
              and returned by \code{\link{ergm.bounddeg}}
              }
        }
      }
 }
}
\seealso{\code{\link{InitMHP}}}
%\examples{}
\keyword{models}