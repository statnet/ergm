% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InitErgmTerm.R, R/InitWtErgmTerm.R
\name{edges-ergmTerm}
\alias{edges-ergmTerm}
\alias{InitErgmTerm.edges}
\alias{InitWtErgmTerm.nonzero}
\alias{nonzero-ergmTerm}
\alias{InitWtErgmTerm.edges}
\title{Number of edges in the network}
\usage{
# binary: edges

# valued: nonzero

# valued: edges
}
\description{
This term adds one network statistic equal to the number of
edges (i.e. nonzero values) in the network. For undirected networks, \code{edges}
is equal to \code{kstar(1)}; for directed networks, edges is equal to both
\code{ostar(1)} and \code{istar(1)}.
}
\seealso{
\code{\link{ergmTerm}} for index of model terms currently visible to the package.

\Sexpr[results=rd,stage=render]{ergm:::.formatTermKeywords("ergmTerm", "edges", "subsection")}
}
\concept{directed}
\concept{dyad-independent}
\concept{undirected}
