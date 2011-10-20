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
\usage{stergm(formation, dissolution, theta.form0 = NULL, theta.diss = NULL, seed = NULL, MH.burnin = 1000, constraints = ~., meanstats = NULL, stergm.order = "FormAndDiss", control = control.stergm(), verbose = FALSE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{formation}{
    An \code{\link{ergm}}-style formula for the formation phase, with
    the LHS network being the observed network or a "template" for it if
    meanstats are used.
}
  \item{dissolution}{
    An \code{\link{ergm}}-style formula for the dissolution phase. The
    LHS is not used.
}
  \item{theta.form0}{
    vector; the starting vector of formation phase parameters for the estimation.
    By default (\code{NULL}), it's estimated using the method specified
    in \code{\link[=control.stergm]{control$initialfit}}.
}
  \item{theta.diss}{
    vector; the dissolution phase parameters. Currently, assumed fixed
    and known.
}
  \item{seed}{
    integer; random number integer seed.  Defaults to
   \code{NULL} to use whatever the state of the random number
   generator is at the time of the call.
}
  \item{MH.burnin}{
    Number of Metropolis-Hastings (one-toggle) steps within each
  formation/dissoution phase. That is, there are \code{2*MH.interval}
  proposals per dynamic time step.
}
  \item{constraints}{
    Constraints for the post-formation and post-dissolution
  networks. See help for \code{\link{ergm}}. \emph{Might not work at
  this time, and have a weird interpretation.}
}
  \item{meanstats}{
    vector; The target equilibrium mean for the network statistics of
  interest (only \code{formation} statistics at this time).
    If it is missing then mean-value parameter used are the observed
    statistics of the network in the formula.
  }
    \item{stergm.order}{
    The "order" in which formation and dissolution phases occur:
    "FormThenDiss" and "DissThenForm" for sequential models, "FormAndDiss"
    for the simultaneous model (the default), and "FormOnly" and
    "DissOnly" for debugging.
  }
  \item{control}{
    A list of control parameters for algorithm
    tuning. Constructed using \code{\link{control.stergm}}. 
  }
  \item{verbose}{
    logical or integer; if positive, the program will print out progress
    information. Higher values result in more output.
  }
  \item{\dots}{
    Additional
    arguments, to be passed to lower-level functions in the future.
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