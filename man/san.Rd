%  File ergm/man/san.Rd
%  Part of the statnet package, http://statnet.org
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
%  Copyright 2012 the statnet development team
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{san}
\alias{san.ergm}
\alias{san.formula}
\alias{san.default}
\alias{san}
\title{Use Simulated Annealing to attempt to match a network to a vector
of mean statistics}
\usage{
\method{san}{formula}(object, 
                      constraints=~., 
                      target.stats=NULL,
                      nsim=1, 
                      basis=NULL,
                      sequential=TRUE,
                      control=control.san(),
                      verbose=FALSE, 
                      \dots)
\method{san}{ergm}(object, 
                   formula=object$formula, 
                   constraints=object$constraints, 
                   target.stats=object$target.stats,
                   nsim=1, 
                   basis=NULL,
                   sequential=TRUE, 
                   control=object$SAN.control,
                   verbose=FALSE, 
                   \dots)
}
\description{
  This function attempts to find a network or networks
  whose statistics match those passed
  in via the \code{target.stats} vector.
}
\arguments{
  \item{object}{Either
    a \code{\link{formula}} or an \code{\link{ergm}} object.
    The \code{\link{formula}} should be of the form \code{y ~ <model terms>},
    where \code{y} is a network object or a matrix that can be
    coerced to a \code{\link[network]{network}} object.  For the details 
    on the possible
    \code{<model terms>}, see \code{\link{ergm-terms}}.  To create a
    \code{\link[network]{network}} object in \R, use the 
    \code{network()} function,
    then add nodal attributes to it using the \code{\%v\%}
    operator if necessary.
  }
  \item{formula}{(By default, the \code{formula}
    is taken from the \code{ergm} object.
    If a different \code{formula} object is wanted, specify it here.
  }
  \item{constraints}{A one-sided formula specifying one or more constraints
    on the support of the distribution of the networks being
    simulated. See the documentation for a similar argument for
    \code{\link{ergm}} for more information. For
    \code{simulate.formula}, defaults to no constraints. For
    \code{simulate.ergm}, defaults to using the same constraints as
    those with which \code{object} was fitted.
  }
  \item{target.stats}{A vector of the same length as the number of terms 
    implied by the formula, which is either \code{object} itself in the 
    case of \code{san.formula} or \code{object$formula} in the case 
    of \code{san.ergm}.
  }
  \item{nsim}{Number of desired networks.
  }
  \item{basis}{If not NULL, a \code{network} object used to start  
    the Markov chain.  If NULL, this is taken to be the network named 
    in the formula.
  }
  \item{sequential}{Logical: If TRUE, the returned draws always use 
    the prior draw as the starting network; if FALSE, 
    they always use the original network.
  }
  \item{control}{A list of control parameters for algorithm
    tuning; see \code{\link{control.san}}.
  }
  \item{verbose}{Logical:  If TRUE, print out more detailed information 
    as the simulation runs.
  }
  \item{\dots}{Further arguments passed to other functions.
  }
}


\value{
  A network or list of networks that hopefully have network statistics close
  to the \code{target.stats} vector.
}
\keyword{models}
