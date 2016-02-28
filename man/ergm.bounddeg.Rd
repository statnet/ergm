\name{ergm.bounddeg}
\alias{ergm.bounddeg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
initializes the parameters to bound degree during sampling
}
\description{
Not normally called directly by user, \code{ergm.bounddeg} initializes the list of parameters used to bound the degree during the Metropolis Hastings sampling process, and issues warnings if the original network doesn't meet the constraints specified by 'bounddeg'. 
}
\usage{
ergm.bounddeg(bounddeg, nw)
}

\arguments{
  \item{bounddeg}{a list of parameters which may contain the following for a network of size n nodes:
  \itemize{
      \item{attribs: an nxp matrix, where entry ij is TRUE if node i has
               attribute j, and FALSE otherwise; default=an nx1 matrix
               of 1's}
      \item{maxout : an nxp matrix, where entry ij is the maximum number of
               out degrees for node i to nodes with attribute j;
               default=an nxp matrix of the value (n-1)}
      \item{maxin  : defined similarly to maxout, but ignored for undirected}
               networks; default=an nxp matrix of the value (n-1)
      \item{minout : defined similarly to maxout; default=an nxp matrix of 0's}
      \item{minin  : defined similarly to maxout, but ignored for undirected
               networks; default=an nxp matrix of 0's}
}
}
  \item{nw}{
the orginal \code{network} specified to \code{ergm} in 'formula'
}
}
\details{
  In some modeling situations, the degree of certain nodes are constrained to
  lie in a certain range (rather than their theoretically possible range
  of 0 to n-1).  Such sample space constraints may be incorporated into the
  ergm modeling process, and if so then the MCMC routine is prevented from 
  visiting network states that violate any of these bounds.  
  
  In case there are categories of nodes and degree bounds for each set of 
  categories, such constraints may be incorporated as well.  For instance,
  if the nodes are girls and boys, and there is a maximum of 5 out-ties to boys and a maximum of 5 out-ties to girls for each node, we would define p to be 2, and the nxp matrix attribs would have TRUE in the first column (say) for exactly those nodes that are boys and TRUE in the second column for only the girls.  The maxout matrix would consist of all 5s in this case, and the other arguments would be left as their default values.
  
  Since the observed network is generally the beginning of the Markov chain, it must satisfy all of the degree constraints itself; thus, this function returns an error message if any bound is violated by the observed network.
}
\value{
  a list of parameters used to bound degree during sampling
  \itemize{
    \item condAllDegExact: always FALSE
    \item attribs        : as defined above
    \item maxout         : as defined above
    \item maxin          : as defined above
    \item minout         : as defined above
    \item minin          : as defined above
}
}


\seealso{
\code{\link{ergm_MH_proposals}}
}
