% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InitErgmTerm.R
\name{b1starmix-ergmTerm}
\alias{b1starmix-ergmTerm}
\alias{InitErgmTerm.b1starmix}
\title{Mixing matrix for \eqn{k}-stars centered on the first mode of a bipartite network}
\usage{
# binary: b1starmix(k, attr, base=NULL, diff=TRUE)
}
\arguments{
\item{k}{only a single value of \eqn{k} is allowed}

\item{attr}{a vertex attribute specification (see Specifying Vertex attributes and Levels (\code{?nodal_attributes}) for details.)}

\item{base}{deprecated}

\item{diff}{whether a different statistic is created for each value seen in a b2 node. When \code{diff=TRUE},
the default, a different statistic is created for each value and thus the behavior of this term is reminiscent of the
\code{nodemix} term, from which it takes its name; when \code{diff=FALSE} ,
all homophilous \eqn{k}-stars are counted together, though these \eqn{k}-stars are still
categorized according to the value of the central b1 node.}
}
\description{
This term counts all \eqn{k}-stars in which
the b2 nodes (called events in some contexts) are homophilous in the sense
that they all share the same value of \code{attr} . However, the b1 node
(in some contexts, the actor) at the center of the \eqn{k}-star does NOT have to
have the same value as the b2 nodes; indeed, the values taken by the b1
nodes may be completely distinct from those of the b2 nodes, which allows
for the use of this term in cases where there are two separate nodal
attributes, one for the b1 nodes and another for the b2 nodes (in this case,
however, these two attributes should be combined to form a single nodal
attribute, \code{attr}). A different statistic is created for each
value of \code{attr} seen in a b1 node, even if no \eqn{k}-stars are observed
with this value.
}
\note{
The argument \code{base} is retained for backwards compatibility and may be
removed in a future version. When both \code{base} and \code{levels} are passed,
\code{levels} overrides \code{base}.
}
\seealso{
\code{\link{ergmTerm}} for index of model terms currently visible to the package.

\Sexpr[results=rd,stage=render]{ergm:::.formatTermKeywords("ergmTerm", "b1starmix", "subsection")}
}
\concept{bipartite}
\concept{categorical nodal attribute}
\concept{undirected}
