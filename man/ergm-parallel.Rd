%  File man/ergm-parallel.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2014 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm-parallel}
\alias{ergm-parallel}
\alias{parallel}
\alias{ergm.parallel}
\alias{parallel.ergm}
\alias{parallel-ergm}

\docType{methods}
\title{Parallel Processing in the \code{\link[=ergm-package]{ergm}} Package}
\description{
  For estimation that require MCMC, \code{\link[=ergm-package]{ergm}} 
  can take advantage of multiple
  CPUs or CPU cores on the system on which it runs, as well as computing
  clusters. It uses package \code{parallel} and \code{snow} to
  facilitate this, and supports all cluster types that they does.

  The number of nodes used and the parallel API are controlled using the
  \code{parallel} and \code{parallel.type} arguments passed to the
  control functions, such as \code{\link{control.ergm}}.
}

% rpvm package is currently not maintained.
%\section{PVM}{
%  To use PVM to accelerate ERGM sampling, pass the control parameter
%  \code{parallel.type="PVM"}.
%  \code{\link[=ergm-package]{ergm}} and \code{\link[snow]{snow}} use
%  \code{\link[rpvm]{rpvm}} package to communicate with a PVM
%  cluster. \code{\link{ergm}} will check if an PVM daemon is already
%  running and will attempt to start one if not.
%}

\section{PSOCK clusters}{
  The \code{parallel} package is used with PSOCK clusters by default,
  to utilize multiple cores on a system. The number of
  cores on a system can be determined with the \code{detectCores} function.
  
  This method works with the base
  installation of R on all platforms, and does not require additional software. 
  
  For more advanced applications, such as clusters that span multiple 
  machines on a network, the clusters can be initialized manually, 
  and passed into \code{ergm} using the
  \code{parallel} control argument. See the second example below.
}

\section{MPI clusters}{
  To use MPI to accelerate ERGM sampling, pass the control parameter
  \code{parallel.type="MPI"}.  
  \code{\link[=ergm-package]{ergm}} requires the \code{snow} and
  \code{Rmpi} packages to communicate with an MPI
  cluster. 
  
  Using MPI clusters requires the system to have an existing MPI installation.
  See the MPI documentation for your particular platform for instructions.

  To use {{{ergm}}} across multiple machines in a high performance computing environment,
  see the section "User initiated clusters" below.
}

\section{User initiated clusters}{
  A cluster can be passed into \code{ergm} with the \code{parallel} control parameter. 
  \code{ergm} will detect the number of nodes in the cluster, and use all of them for
  MCMC sampling. This method is flexible: it will accept any cluster type that is
  compatible with \code{snow} or \code{parallel} packages. Usage examples
  for a multiple-machine high performance MPI cluster can be found at the statnet wiki: 
\url{https://statnet.csde.washington.edu/trac/wiki/ergmParallel} 
}

\examples{
\donttest{
# Uses 2 SOCK clusters for MCMLE estimation
data(faux.mesa.high)
nw <- faux.mesa.high
fauxmodel.01 <- ergm(nw ~ edges + isolates + gwesp(0.2, fixed=TRUE), 
                     control=control.ergm(parallel=2, parallel.type="PSOCK"))
summary(fauxmodel.01)

}
}
