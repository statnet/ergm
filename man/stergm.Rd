%  File ergm/man/stergm.Rd
%  Part of the statnet package, http://statnet.org
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
%  Copyright 2012 the statnet development team
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{stergm}
\alias{stergm}
\alias{stergm.object}
\title{Separable Temporal Exponential Family Random Graph Models
}
\description{
  \code{\link{stergm}} is used for finding Separable Temporal ERGMs'
  (STERGMs) Conditional MLE (CMLE) (Krivitsky and Handcock, 2010) and
  Equilibrium Generalized Method of Moments Estimator (EGMME) (Krivitsky, 2009).
}
\usage{
  stergm (nw, 
          formation, 
          dissolution, 
          estimate, 
          times=NULL, 
          offset.coef.form=NULL, 
          offset.coef.diss=NULL,
          targets=NULL, 
          target.stats=NULL,
          eval.loglik=FALSE,
          control=control.stergm(),
          verbose=FALSE, 
          \dots)
}
\arguments{
  \item{nw}{A \code{\link[network]{network}} object (for EGMME); or
    \code{\link[networkDynamic]{networkDynamic}} object, a
    \code{\link{network.list}} object, or a \code{\link{list}}
    containing networks (for CMLE and CMPLE).
  }
  \item{formation, dissolution}{
    One-sided \code{\link{ergm}}-style formulas for the formation and
    dissolution models, respectively.
  }
  
  \item{estimate}{One of "EGMME" for Equilibrium Generalized Method
    of Moments Estimation, based on a single network with some temporal
    information and making an assumption that it is a product of a
    STERGM process running to its stationary (equilibrium) distribution; "CMLE" 
    for Conditional Maximum Likelihood Estimation, modeling a transition
    between two networks, or
    "CMPLE" for Conditional Maximum PseudoLikelihood Estimation, using
    MPLE instead of MLE.
    CMPLE is extremely inaccurate at this time.
  }
  \item{times}{For CMLE and CMPLE estimation, times or indexes at which
    the networks whose transition is to be modeled are
    observed. Default to \code{c(0,1)} if \code{nw} is a
    \code{\link[networkDynamic]{networkDynamic}} and to \code{c(1:2)}
    (first two elements) if \code{nw} is a \code{\link{network.list}} or
    a \code{\link{list}}. Unused for EGMME.
  }
  \item{offset.coef.form}{Numeric vector to specify offset formation
    parameters.
  }
  \item{offset.coef.diss}{Numeric vector to specify offset dissolution
    parameters.
  }
  \item{targets}{One-sided \code{\link{ergm}}-style formula specifying statistics whose
    moments are used for the EGMME. Unused for CMLE and CMPLE.
  }

  \item{target.stats}{
    A vector specifying the values of the \code{targets} statistics that
    EGMME will try to match. 
    Defaults to the
    statistics of \code{nw}. Unused for CMLE and CMPLE.
  }
  \item{eval.loglik}{Whether or not to calculate the log-likelihood of
    a CMLE STERGM fit. See \code{\link{ergm}} for details.}
  \item{control}{
    A list of control parameters for algorithm
    tuning. Constructed using \code{\link{control.stergm}}. 
  }
  \item{verbose}{
    logical or integer; if TRUE or positive, the program will print out progress
    information. Higher values result in more output.
  }
  \item{\dots}{Additional arguments, to be passed to lower-level 
    functions.
  }  
}
\section{Model Terms}{
  See \code{\link{ergm}} and  \code{\link{ergm-terms}} for details. At
  this time, only linear ERGM terms are allowed.
}
\value{\code{\link{stergm}} returns an object of class \code{\link{stergm}} that is a list
  consisting of the following elements:
  \item{formation, dissolution}{Formation and dissolution formulas,
    respectively.}
  \item{targets}{The targets formula.}
  \item{target.stats}{The target statistics.}
  \item{estimate}{The type of estimate.}
  \item{opt.history}{A matrix containing the full trace of the EGMME
    optimization process: coefficients tried and target statistics
    simulated.}
  \item{sample}{An \code{\link{mcmc}} object containing target
    statistics sampled at the estimate.}
  \item{covar}{The full estimated variance-covariance matrix of the
    parameter estimates for EGMME. (Note that although the CMLE
    formation parameter estimates are independent of the dissolution
    parameter estimates due to the separability assumption, this is not
    necessarily the case for EGMME.) }
  \item{formation.fit, dissolution.fit}{For CMLE and CMPLE,
    \code{\link{ergm}} objects from fitting formation and dissolution,
    respectively. For EGMME, stripped down \code{\link{ergm}}-like lists.}
  \item{network}{The original network.}
  \item{control}{The control parameters used to fit the model.}
  
  See the method \code{\link{print.stergm}} for details on how
  an \code{\link{stergm}} object is printed.  Note that the
  method \code{\link{summary.stergm}} returns a summary of the
  relevant parts of the \code{\link{stergm}} object in concise summary
  format.
}

\references{
  Krivitsky PN (2009). Statistical Models for Social Network Data and
  Processes. PhD Thesis. University of Washington.

  Krivitsky PN, Handcock MS (2010). A Separable Model for Dynamic
  Networks. \url{http://arxiv.org/abs/1011.1937}
}

\seealso{network, \%v\%, \%n\%, \code{\link{ergm-terms}}}