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
\method{logLik}{ergm}(object, nsteps = 20, add = FALSE, eval.loglik = add, force.reeval = FALSE, ...)
}
\arguments{
  \item{object}{
An \code{\link[=ergm.object]{ergm}} fit, returned by \code{\link{ergm}}.
}
  \item{nsteps}{
Number of bridges for bridge sampling estimation.
}
  \item{add}{
Instead of returning the log-likelihood, return \code{object} with
log-likelihood value set.
}
  \item{eval.loglik}{
Evaluate the log-likelihood if not set on \code{object}.
}
  \item{force.reeval}{
Whether to reestimate the log-likelihood even if \code{object} already
has an estiamte. Defaults to \code{FALSE}.
}
  \item{\dots}{
Other arguments to the bridge sampling functions.
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
\author{
Pavel N. Krivitsky
}

\seealso{
\code{\link{logLik}}, \code{\link{ergm.bridge.llr}}, \code{\link{ergm.bridge.dindstart.llk}}
}
\examples{
# See help(ergm) for a description of this model:
data(florentine)
gest <- ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle)
# Log-likelihood is not evaluated, so no deviance, AIC, or BIC:
summary(gest)
# Evaluate the log-likelihood and attach it to the object.
gest <- logLik(gest, add=TRUE)
# Deviances, AIC, and BIC now shown:
summary(gest)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models  }
