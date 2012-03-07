%  File ergm/man/logLik.ergm.Rd
%  Part of the statnet package, http://statnet.org
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
%  Copyright 2012 the statnet development team
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{logLik.ergm}
\alias{logLik.ergm}
\title{
  A \code{\link{logLik}} method for \code{\link[=ergm.object]{ergm}}.
}
\description{
A function to return the log-likelihood associated with an
\code{\link[=ergm.object]{ergm}} fit, evaluating it if necessary.
}
\usage{
\method{logLik}{ergm}(object, 
                      add=FALSE, 
                      force.reeval=FALSE, 
                      eval.loglik=add || force.reeval, 
                      control=control.logLik.ergm(), 
                      \dots)
}
\arguments{
  \item{object}{
    An \code{\link[=ergm.object]{ergm}} fit, returned by \code{\link{ergm}}.
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
  If \code{add=FALSE} (the default), a \code{\link{logLik}} object. If
  \code{add=TRUE} (the default), an \code{\link[=ergm.object]{ergm}}
  object with the log-likelihood set.
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
# See help(ergm) for a description of this model:
data(florentine)
gest <- ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle, eval.loglik=FALSE)
# Log-likelihood is not evaluated, so no deviance, AIC, or BIC:
summary(gest)
# Evaluate the log-likelihood and attach it to the object.
gest <- logLik(gest, add=TRUE)
# Deviances, AIC, and BIC now shown:
summary(gest)
}
\keyword{ models  }
