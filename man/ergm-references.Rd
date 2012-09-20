\name{ergm-references}
\alias{ergm-references}
\alias{references-ergm}
\alias{ergm.references}
\alias{references.ergm}
\alias{Bernoulli}
\alias{InitReference.Bernoulli}

\docType{methods}

\title{Reference Measures for Exponential-Family Random Graph Models
  (binary ERGMs)}

\description{
  This page documents the default (Bernoulli) reference measure for
  binary ERGMs. It does not need to be specified explicitly, so it is
  included mainly for completeness.

  The reference measure is specified on the RHS of a one-sided formula passed as
  the \code{reference} argument to \code{\link[ergm]{ergm}}. 
  See the \code{\link[ergm]{ergm}} documentation for a complete
  description of how reference measures are specified.
}

\section{Possible reference measures to represent binary ERGMs}{

  Reference measures currently available are:
  \describe{
    
    \item{\code{Bernoulli}}{\emph{Bernoulli-reference ERGM:} Specifies each
      dyad's baseline distribution to be Bernoulli with probability of
      the tie being \eqn{0.5}. }
    }
}

\references{

Hunter DR, Handcock MS, Butts CT, Goodreau SM, Morris M (2008b).
\pkg{ergm}: A Package to Fit, Simulate and Diagnose
Exponential-Family Models for Networks.
\emph{Journal of Statistical Software}, 24(3).
\url{http://www.jstatsoft.org/v24/i03/}.
  
Krivitsky PN (2012). Exponential-Family Random Graph Models for Valued Networks. \emph{Electronic Journal of Statistics}, 2012, 6, 1100-1128. \href{http://dx.doi.org/10.1214/12-EJS696}{\code{doi:10.1214/12-EJS696}}
}

\seealso{ergm, network, \%v\%, \%n\%, sna,
         summary.ergm, print.ergm}
\keyword{models}