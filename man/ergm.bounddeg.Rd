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
%%  ~~ If necessary, more details than the description above ~~
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
