%  File man/ergm-parallel.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2013 Statnet Commons
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
  \code{\link[=ergm-package]{ergm}} can take advantage of multiple
  CPUs or CPU cores on the system on which it runs, as well as computing
  clusters. It uses package \code{snow} to
  facilitate this, and supports all cluster types that it does.

  The parallel API and the number of nodes used are controlled using the
  \code{parallel} and \code{parallel.type} arguments passed to the
  control functions, such as \code{\link{control.ergm}}.

  This entry describes common problems and workarounds associated with
  particular parallel processing APIs.  }

% rpvm package is currently not maintained.
%\section{PVM}{
%  To use PVM to accelerate ERGM sampling, pass the control parameter
%  \code{parallel.type="PVM"}.
%  \code{\link[=ergm-package]{ergm}} and \code{\link[snow]{snow}} use
%  \code{\link[rpvm]{rpvm}} package to communicate with a PVM
%  cluster. \code{\link{ergm}} will check if an PVM daemon is already
%  running and will attempt to start one if not.
%}

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
  type="MPI")}). \code{\link{ergm}} will notice the
  preexisting cluster and make use of it, but it will not stop it.

}

\examples{
\donttest{
# See help(ergm) for a description of this model.
data(florentine)
gest <- ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangles,
             eval.loglik=FALSE,
             control=control.ergm(parallel=2, parallel.type="SOCK"))
summary(gest)
# Note the combined MCMC diagnostics:
mcmc.diagnostics(gest)
}
}
