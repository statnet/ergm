\name{check.ErgmTerm}
\alias{check.ErgmTerm}
\alias{ergm.checkargs}
\alias{ergm.checkbipartite}
\alias{ergm.checkconstraints.model}
\alias{ergm.checkdegeneracy}
\alias{ergm.checkdirected}
\alias{ergm.checkextreme.model}
\title{
Ensures an Ergm Term and its Arguments Meet Appropriate Conditions
}
\description{
This is a low-level function not intended to be called by end users.  For information on ergm terms, see \link{ergm-terms} The \code{check.ErgmTerm} function ensures for the \code{InitErgmTerm.X} function that the term X:
\itemize{
  \item is applicable given the 'directed' and 'bipartite' attributes of the given network
  \item is not applied to a directed bipartite network
  \item has an appropiate number of arguments
  \item has correct argument types if arguments where provided
  \item has default values assigned if defaults are available
   }
 by halting execution if any of the first 3 criteria are not met.
}
\usage{
check.ErgmTerm(nw, arglist, directed = NULL, bipartite = NULL,
               nonnegative = FALSE, varnames = NULL, vartypes = NULL, 
               defaultvalues = list(), required = NULL, response = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nw}{
the network that term X is being checked against  
}
  \item{arglist}{
the list of arguments for term X
}
  \item{directed}{
whether term X requires a directed network (T or F); default=NULL
}
  \item{bipartite}{
whether term X requires a bipartite network (T or F); default=NULL
}
  \item{nonnegative}{
whether term X requires a network with only nonnegative weights; default=FALSE
}
  \item{varnames}{
the vector of names of the possible arguments for term X; default=NULL 
}
  \item{vartypes}{
the vector of types of the possible arguments for term X; default=NULL 
}
  \item{defaultvalues}{
the list of default values for the possible arguments of term X; default=list()
}
  \item{required}{
the logical vector of whether each possible argument is required;  default=NULL
}
  \item{response}{
edge attribute name?
}
}

\value{
a list of the values for each possible argument of term X; user provided  values are used when given, default values otherwise.
}

