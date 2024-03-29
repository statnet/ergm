% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InitErgmTerm.R
\name{b2star-ergmTerm}
\alias{b2star-ergmTerm}
\alias{InitErgmTerm.b2star}
\title{\eqn{k}-stars for the second mode in a bipartite network}
\usage{
# binary: b2star(k, attr=NULL, levels=NULL)
}
\arguments{
\item{k}{a vector of distinct integers}

\item{attr, levels}{a vertex attribute specification; if \code{attr} is
specified, then the count is over the instances where all nodes
involved have the same value of the attribute. \code{levels} specified
which values of \code{attr} are included in the count. (See Specifying
Vertex attributes and Levels
(\code{\link[=nodal_attributes]{?nodal_attributes}}) for details.)}
}
\description{
This term adds one network statistic to the model for
each element in \code{k} . The \eqn{i} th such statistic counts the number of
distinct \code{k[i]} -stars whose center node is in the second mode of the
network. The second mode of a bipartite network object is sometimes known as
the "event" mode. A \eqn{k} -star is defined to be a center node \eqn{N} and
a set of \eqn{k} different nodes \eqn{\{O_1, \dots, O_k\}}{\{O[1], ..., O[k]\}} such that the
ties \eqn{\{N, O_i\}} exist for \eqn{i=1, \dots, k} . This term can only be used for
undirected bipartite networks.
}
\note{
\code{b2star(1)} is equal to \code{b1star(1)} and to \code{edges} .
}
\seealso{
\code{\link{ergmTerm}} for index of model terms currently visible to the package.

\Sexpr[results=rd,stage=render]{ergm:::.formatTermKeywords("ergmTerm", "b2star", "subsection")}
}
\concept{bipartite}
\concept{categorical nodal attribute}
\concept{undirected}
