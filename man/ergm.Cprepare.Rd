\name{ergm.Cprepare}
\alias{ergm.Cprepare}
\alias{ergm.Cprepare.el}
\alias{ergm.Cprepare.miss}
\alias{ergm.design}

\title{
Internal Function to Prepare Data for ergm's C Interface 
}
\description{
These are internal functions not intended to be called by end users.  The \code{ergm.Cprepare} function builds an object called Clist that contains all the necessary ingredients to be passed to the C functions, other functions create edgelists and handle missing edge data. 
}
\usage{
ergm.Cprepare(nw, m, response = NULL)

ergm.Cprepare.el(x, attrname=NULL, prototype=NULL)

ergm.Cprepare.miss(nw)

ergm.design(nw, model, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nw,x}{
a network object
}
  \item{m,model}{
a model object, as returned by \code{\link{ergm.getmodel}}
}
  \item{response,attrname}{
character name of an edge attribute
} 
\item{prototype}{
  A network whose relevant attributes (size, directedness, bipartitedness,
  and presence of loops) are imposed on the output edgelist if \code{x}
  is already an edgelist. (For example, if the prototype is undirected,
  \code{ergm.Cprepare.el} will ensure that \eqn{t < h}.)
}
\item{verbose}{logical, whether the design matrix should be printed; default=FALSE}
}
\details{
These low-level functions are used by other ergm-related packages, but should never need to be called directly by the user. 
\itemize{
\item \code{ergm.Cprepare} builds an object called Clist that contains all the necessary ingredients to be passed to the C functions
\item \code{ergm.Cprepare.el} constructs and serializes a very simple static edgelist, with the vertex having the lesser index the tail and sorted by tails, then by heads.
\item \code{ergm.Cprepare.miss} constructs an edgelist as \code{ergm.Cprepare.el}, but only includes 'missing' edges (marked as NA)
\item \code{ergm.design}  functions as \code{ergm.Cprepare} would, but acts on the network of missing edges
}
}
\value{
\code{ergm.Cprepare} returns Clist:  a list of parameters used by several of the fitting routines containing
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
 
\code{ergm.design} returns: Clist.miss
\itemize{
\item{if 'nw' has missing edges, see the return list, 'Clist', from the
                              \code{ergm.Cprepare} function header}
\item{if 'nw' hasn't any missing edges, the list will only contain NULL
                                values for the 'tails' and 'heads' components,
                                a 0 for 'nedges' and 'dir' appropriately set }
}
\code{ergm.Cprepare.miss} returns a vector of length 1+Nmissing*2. The first element is the number of missing edges, and the remainder a column-major edgelist
}
