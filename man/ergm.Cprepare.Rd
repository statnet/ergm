\name{ergm.Cprepare}
\alias{ergm.Cprepare}
\alias{ergm.Cprepare.el}
\alias{ergm.Cprepare.miss}

\title{
Internal Function to Prepare Data for ergm's C Interface 
}
\description{
This is an internal function not intended to be called by end users.  The \code{ergm.Cprepare} function builds an object called Clist that contains all the necessary ingredients to be passed to the C functions
}
\usage{
ergm.Cprepare(nw, m, response = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nw}{
a network object
}
  \item{m}{
a model object, as returned by \code{\link{ergm.getmodel}}
}
  \item{response}{
character name of an edge attribute
}
}
\details{
This function is used directly by other ergm-related packages, but should never need to be called directly by the user. 
}
\value{
Clist:  a list of parameters used by several of the fitting routines containing
\itemize{
\item{      n           :  the size of the network}
\item{     dir         :  whether the network is directed (T or F)}
\item{     bipartite   :  whether the network is bipartite (T or F)}
\item{     ndyads      :  the number of dyads in the network}
\item{     nedges      :  the number of edges in this network}
\item{     tails       :  the vector of tail nodes; tail nodes are
                               the 1st column of the implicit edgelist,
                               so either the lower-numbered nodes in an
                               undirected graph, or the out nodes of a
                               directed graph, or the b1 nodes of a bi-
                               partite graph}
\item{     heads           :  the vector of head nodes; head nodes are
                               the 2nd column of the implicit edgelist,
                               so either the higher-numbered nodes in an
                               undirected graph, or the in nodes of a
                               directed graph, or the b2 nodes of a bi-
                               partite graph}
\item{     nterms      :  the number of model terms}
\item{     nstats      :  the total number of change statistics
                           for all model terms}
\item{     inputs      :  the concatenated vector of 'input's from each
                           model term as returned by <InitErgmTerm.X> or
                           <InitErgm.X>}
\item{     fnamestring :  the concatenated string of model term names}
\item{     snamestring :  the concatenated string of package names that
                           contain the C function 'd_fname'; default="ergm"
                           for each fname in fnamestring}
 }
}
