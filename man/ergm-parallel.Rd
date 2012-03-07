%  File ergm/man/ergm-parallel.Rd
%  Part of the statnet package, http://statnet.org
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
%  Copyright 2012 the statnet development team
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm-parallel}
\alias{ergm-parallel}
\alias{parallel}
\alias{ergm.parallel}
\alias{parallel.ergm}
\alias{parallel-ergm}

\docType{methods}
\title{Parallel Processing in the \code{\link[=ergm-package]{ergm}} Package}
\description{
  \code{\link[=ergm-package]{ergm}} can take advantage of multiple
  CPUs or CPU cores on the system on which it runs, as well as computing
  clusters.

  This entry describes the APIs it can use, as well as common problems
  and workarounds associated with parallel processing in \code{\link[=ergm-package]{ergm}}.
}


\section{MPI}{
  To use MPI to accellerate ERGM sampling, pass the control parameter
  \code{parallel.type="MPI"}.
  \code{\link[=ergm-package]{ergm}} and \code{snow} use
  \code{Rmpi} package to communicate with an MPI
  cluster. \code{\link{ergm}} will check if an MPI cluster already
  exists and will create one if one doesn't.

  On some installations, the function \code{\link[snow]{stopCluster}}
  does not work properly for MPI clusters. Because \code{\link{ergm}}
  creates a cluster every time it needs an MCMC sample and disbands it
  once it finishes, using MPI on these installations may fail. A
  workaround is to start the cluster outside of \code{\link{ergm}}
  (e.g., \code{dummy <- makeCluster(nnodes,
  type="MPI")}). \code{\link{ergm}} will notice the preexisting
  preexisting cluster and make use of it, but it will not stop it.

}