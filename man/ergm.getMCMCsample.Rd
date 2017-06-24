%  File man/ergm.getMCMCsample.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2017 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm.getMCMCsample}
\alias{ergm.getMCMCsample}
\alias{ergm.mcmcslave}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal Function to Sample Networks Using C Wrapper
}
\description{
This is an internal function, not normally called directly by the user. The \code{ergm.getMCMCsample} function samples networks using an MCMC algorithm via
\code{MCMC_wrapper.C} and is caple of running in multiple threads using \code{ergm.mcmcslave}.  

 The \code{ergm.mcmcslave} function is that which the slave nodes in a parallel process will call to perform  a validation on the mcmc equal to their slave number. It also returns an MCMC sample.


}
\usage{
ergm.getMCMCsample(nw, model, MHproposal, eta0, control, verbose, response = NULL, ...)

ergm.mcmcslave(Clist, MHproposal, eta0, control, verbose, ..., prev.run = NULL, 
    burnin = NULL, samplesize = NULL, interval = NULL, maxedges = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nw}{
a network object
}
  \item{model}{
a model for the given 'nw' as returned by <ergm.getmodel>
}
  \item{MHproposal}{
a list of the parameters needed for Metropolis-Hastings proposals and
                the result of calling <MHproposal>
}
  \item{eta0}{
 the initial eta coefficients
}
  \item{control}{
list of MCMC tuning parameters; (see \link{control.ergm}) 
}
  \item{verbose}{
whether the C functions should be verbose; default=FALSE
}
  \item{response}{
characher, name of an edge attribute
}
  \item{\dots}{
additional arugments?
}

\item{Clist}{the list of parameters returned by \code{\link{ergm.Cprepare}}}
\item{prev.run}{ output stats from previous run ??? }
\item{burnin}{ number of proposals before any MCMC sampling is done. see \link{control.ergm} MCMC.burnin}
\item{samplesize}{number of network statistics, randomly drawn from a given distribution on the set of all networks, returned by the Metropolis-Hastings algorithm. see \link{control.ergm} MCMC.samplesize}
\item{interval}{number of proposals between sampled statistics. see \link{control.ergm} MCMC.interval} 
\item{maxedges}{maximum number of edges expected in network. see \link{control.ergm} MCMC.init.maxedges}

}
\details{
Note that the returned stats will be relative to the original network, i.e., the calling function must shift the statistics if required. The calling function must also attach column names to the statistics matrix if required.
}
\value{
   for \code{ergm.getMCMCsample},the sample as a list containing:
   \itemize{
     \item{statsmatrix:  the stats matrix for the sampled networks, RELATIVE TO THE ORIGINAL
                   NETWORK!}
     \item{newnetwork :  the edgelist of the final sampled network}
     \item{nedges     :  the number of edges in the 'newnetwork'} 
     }
     for \code{ergm.mcmcslave} the MCMC sample as a list of the following:
        \itemize{
\item{    s         : the statsmatrix}
\item{     newnwtails: the vector of tails for the new network- is this the final
                 network sampled? - is this the original nw if 'maxedges' is 0}
\item{     newnwheads: the vector of heads for the new network - same q's}
}
     
}

