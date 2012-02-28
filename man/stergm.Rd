\name{stergm}
\alias{stergm}
\alias{stergm.object}
\title{Separable Transition Exponential Family Random Graph Models
}
\description{
  \code{\link{stergm}} is used for finding equilibrium GMoME of Separable Transition ERGMs,
  described by Krivitsky (2009)  A number of methods is possible, with a
  very limited Robbins-Monro being the default.
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
  \item{nw}{A \code{\link[network]{network}} object giving the
    properties of the sample space of networks, such as the size,
    directedness, and covariates.
  }
  \item{formation}{
    An \code{\link{ergm}}-style formula for the formation phase, with
    the left-hand-side network being the observed network or a 
    "template" for it if target.stats are used.
  }
  \item{dissolution}{
    An \code{\link{ergm}}-style formula for the dissolution phase. The
    left-hand-side is not used.
  }
  \item{estimate}{Not yet documented.
  }
  \item{times}{Not yet documented.
  }
  \item{offset.coef.form}{Not yet documented.
  }
  \item{offset.coef.diss}{Not yet documented.
  }
  \item{targets}{Not yet documented.
  }
%  \item{theta.form0}{
%    vector; the starting vector of formation phase parameters for the estimation.
%    By default (\code{NULL}), it's estimated using the method specified
%    in \code{\link[=control.stergm]{control$initialfit}}.
%  }
%  \item{theta.diss}{
%    vector; the dissolution phase parameters. Currently, assumed fixed
%    and known.
%  }
  \item{target.stats}{
    vector; The target equilibrium mean for the network statistics of
    interest (only \code{formation} statistics at this time).
    If it is missing then mean-value parameter used are the observed
    statistics of the network in the formula.
  }
  \item{eval.loglik}{Not yet documented.
  }
  \item{control}{
    A list of control parameters for algorithm
    tuning. Constructed using \code{\link{control.stergm}}. 
  }
  \item{verbose}{
    logical or integer; if positive, the program will print out progress
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
  \item{coef.form, coef.diss}{GMoMEs for formation and dissolution.}
  \item{MH.burnin}{Number of Metropolis-Hastings step during each phase.}
  \item{network}{Original network}
  \item{theta.original,theta.form.original}{The first value of theta.}
 \item{etamap}{The set of functions mapping the true parameter theta
 to the canonical parameter eta (irrelevant except in a curved exponential
 family model)}
 \item{formation, dissolution}{The original \code{\link{formula}}s entered into the \code{\link{stergm}} function.}
 \item{constraints}{Constraints used by original \code{ergm} call}
 \item{prop.weights.form, prop.weights.diss}{MCMC proposal weights used by original \code{stergm} call
 (part of the \code{\link{control.stergm}} function output).}
 \item{offset}{vector of logical telling which model parameters are to be set
 at a fixed value (i.e., not estimated).}
  See the method \code{\link{print.stergm}} for details on how
  an \code{\link{stergm}} object is printed.  Note that the
  method \code{\link{summary.stergm}} returns a summary of the
  relevant parts of the \code{\link{stergm}} object in concise summary
  format.
}

\references{
  Krivitsky PN (2009). Statistical Models for Social Network Data and
  Processes. PhD Thesis. University of Washington.
}

\seealso{network, \%v\%, \%n\%, \code{\link{ergm-terms}}}