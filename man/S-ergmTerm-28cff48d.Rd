% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InitErgmTerm.operator.R,
%   R/InitWtErgmTerm.operator.R
\name{S-ergmTerm}
\alias{S-ergmTerm}
\alias{InitErgmTerm.S}
\alias{InitWtErgmTerm.S}
\title{Evaluation on an induced subgraph}
\usage{
# binary: S(formula, attrs)

# valued: S(formula, attrs)
}
\arguments{
\item{formula}{a one-sided \code{\link[=ergm]{ergm()}}-style formula with the terms to be evaluated}

\item{attrs}{a two-sided formula to be used. A one-sided formula (e.g., \code{~A} ) is symmetrized (e.g., \code{A~A} ).}
}
\description{
This operator takes a two-sided forumla \code{attrs} whose LHS gives the attribute or attribute function for which tails and heads will be used to construct the induced subgraph. They must evaluate either to a logical vector equal in length to the number of tails (for LHS) and heads (for RHS) indicating which nodes are to be used to induce the subgraph or a numeric vector giving their indices.
}
\details{
As with indexing vectors, the logical vector will be recycled to the size of the network or the size of the appropriate bipartition, and negative indices will deselect vertices.

When the two sets are identical, the induced subgraph retains the directedness of the original graph. Otherwise, an undirected bipartite graph is induced.
}
\seealso{
\code{\link{ergmTerm}} for index of model terms currently visible to the package.

\Sexpr[results=rd,stage=render]{ergm:::.formatTermKeywords("ergmTerm", "S", "subsection")}
}
\concept{operator}
