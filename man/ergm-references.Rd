%  File man/ergm-references.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, https://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  https://statnet.org/attribution
%
%  Copyright 2003-2019 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm-references}
\alias{ergm-references}
\alias{references-ergm}
\alias{ergm.references}
\alias{references.ergm}
\alias{Bernoulli}
\alias{Unif}
\alias{DiscUnif}
\alias{StdNormal}
%% \alias{InitErgmReference.Bernoulli}
%% \alias{InitErgmReference.DiscUnif}
%% \alias{InitErgmReference.Unif}
%% \alias{InitErgmReference.StdNormal}
%% \alias{InitWtErgmProposal.DiscUnif}
%% \alias{InitWtErgmProposal.DiscUnifNonObserved}
%% \alias{InitWtErgmProposal.Unif}
%% \alias{InitWtErgmProposal.UnifNonObserved}
%% \alias{InitWtErgmProposal.StdNormal}

\docType{methods}

\title{Reference Measures for Exponential-Family Random Graph Models}

\description{
  This page describes the possible reference measures (baseline distributions)
  for found in the \code{\link[=ergm-package]{ergm}} package, particularly the default (Bernoulli) reference measure for
  binary ERGMs.
  
  The reference measure is specified on the RHS of a one-sided formula passed as
  the \code{reference} argument to \code{\link[ergm]{ergm}}. 
  See the \code{\link[ergm]{ergm}} documentation for a complete
  description of how reference measures are specified.
}

\section{Possible reference measures to represent baseline distributions}{

  Reference measures currently available are:
  \describe{
    
    \item{\code{Bernoulli}}{\emph{Bernoulli-reference ERGM:} Specifies each
      dyad's baseline distribution to be Bernoulli with probability of
      the tie being \eqn{0.5}. This is the only reference measure used
      in binary mode. }

    \item{\code{DiscUnif(a,b)}}{\emph{Discrete-Uniform-reference ERGM:}
      Specifies each dyad's baseline distribution to be discrete uniform
      between \code{a} and \code{b} (both inclusive): \eqn{h(y)=1}, with
      the support being
      \code{a},\code{a}+1,\ldots,\code{b}-1,\code{b}. At this time, both \code{a} and
      \code{b} must be finite. }

    \item{\code{Unif(a,b)}}{\emph{Coninuous-Uniform-reference ERGM:}
      Specifies each dyad's baseline distribution to be continuous uniform
      between \code{a} and \code{b} (both inclusive): \eqn{h(y)=1}, with
      the support being
      [\code{a},\code{b}]. At this time, both \code{a} and
      \code{b} must be finite.  }
      
    \item{\code{StdNormal}}{\emph{Standard-Normal-reference ERGM:}
      Specifies each dyad's baseline distribution to be the normal distribution
      with mean 0 and variance 1.  }

  }
}

\references{

Hunter DR, Handcock MS, Butts CT, Goodreau SM, Morris M (2008b).
\pkg{ergm}: A Package to Fit, Simulate and Diagnose
Exponential-Family Models for Networks.
\emph{Journal of Statistical Software}, 24(3).
\url{https://www.jstatsoft.org/v24/i03/}.
  
Krivitsky PN (2012). Exponential-Family Random Graph Models for Valued Networks. \emph{Electronic Journal of Statistics}, 2012, 6, 1100-1128. \doi{10.1214/12-EJS696}
}

\seealso{ergm, network, \%v\%, \%n\%, sna,
         summary.ergm, print.ergm}
\keyword{models}