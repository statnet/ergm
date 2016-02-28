%  File man/logLik.ergm.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2015 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{logLik.ergm}
\alias{logLik.ergm}
\alias{logLikNull}
\alias{logLikNull.ergm}
\title{
  A \code{\link{logLik}} method for \code{\link[=ergm.object]{ergm}}.
}
\description{
  A function to return the log-likelihood associated with an
  \code{\link[=ergm.object]{ergm}} fit, evaluating it if
  necessary. \code{logLikNull} computes, when possible (see Value), the
  log-probability of observing the observed, unconstrained dyads of the
  network observed under the null model.
}
\usage{
\method{logLik}{ergm}(object, 
                      add=FALSE, 
                      force.reeval=FALSE, 
                      eval.loglik=add || force.reeval, 
                      control=control.logLik.ergm(), 
                      \dots)

logLikNull(object, \dots)

\method{logLikNull}{ergm}(object, 
                          control=control.logLik.ergm(), 
                          \dots)

}
\arguments{
  \item{object}{
    An \code{\link[=ergm.object]{ergm}} fit, returned by
                      \code{\link{ergm}}.
  }
  \item{add}{Logical:  If TRUE, instead of 
    returning the log-likelihood, return \code{object} with
    log-likelihood value set.
  }
  \item{force.reeval}{Logical:  If TRUE,
    reestimate the log-likelihood even if \code{object} already
    has an estiamte. 
  }
  \item{eval.loglik}{Logical:  If TRUE,
    evaluate the log-likelihood if not set on \code{object}.
  }
  \item{control}{
    A list of control parameters for algorithm
    tuning. Constructed using \code{\link{control.logLik.ergm}}. 
  }
  \item{\dots}{
    Other arguments to the likelihood functions.
  }
}
\details{If the log-likelihood was not computed for \code{object},
  produces an error unless \code{eval.loglik=TRUE}
}
\value{
  The form of the output of \code{logLik.ergm} depends on \code{add}: \code{add=FALSE} (the default), a \code{\link{logLik}} object. If
  \code{add=TRUE} (the default), an \code{\link[=ergm.object]{ergm}}
  object with the log-likelihood set.

  \code{logLikNull} returns an object of type
  \code{\link{logLik}} if it is able to compute the null model
  probability, and \code{NA} otherwise.

  As of version 3.1, all likelihoods for which \code{logLikNull} is not
  implemented are computed relative to the reference measure. (I.e., a
  null model, with no terms, is defined to have likelihood of 0, and all
  other models are defined relative to that.)
}
\references{
  Hunter, D. R. and Handcock, M. S. (2006)
  \emph{Inference in curved exponential family models for networks},
  Journal of Computational and Graphical Statistics.
}

\seealso{
  \code{\link{logLik}}, \code{\link{ergm.bridge.llr}}, \code{\link{ergm.bridge.dindstart.llk}}
}
\examples{
# See help(ergm) for a description of this model. The likelihood will
# not be evaluated.
data(florentine)
\dontrun{
# The default maximum number of iterations is currently 20. We'll only
# use 2 here for speed's sake.
gest <- ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle, eval.loglik=FALSE)

gest <- ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle, eval.loglik=FALSE,
             control=control.ergm(MCMLE.maxit=2))
# Log-likelihood is not evaluated, so no deviance, AIC, or BIC:
summary(gest)
# Evaluate the log-likelihood and attach it to the object.

# The default number of bridges is currently 20. We'll only use 3 here
# for speed's sake.
gest.logLik <- logLik(gest, add=TRUE)

gest.logLik <- logLik(gest, add=TRUE, control=control.logLik.ergm(nsteps=3))
# Deviances, AIC, and BIC are now shown:
summary(gest.logLik)
# Null model likelihood can also be evaluated, but not for all constraints:
logLikNull(gest) # == network.dyadcount(flomarriage)*log(1/2)
}
}
\keyword{models}
