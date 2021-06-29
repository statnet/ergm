%  File man/ergm-references.Rd in package ergm, part of the
%  Statnet suite of packages for network analysis, https://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  https://statnet.org/attribution .
%
%  Copyright 2003-2021 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm-references}
\alias{ergm-references}
\alias{references-ergm}
\alias{ergm.references}
\alias{references.ergm}
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
  Reference measures currently available can be accessed via \code{\link{ergmReference}}
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
