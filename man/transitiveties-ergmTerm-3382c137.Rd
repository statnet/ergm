% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InitErgmTerm.transitiveties.R
\name{transitiveties-ergmTerm}
\alias{transitiveties-ergmTerm}
\alias{InitErgmTerm.transitiveties}
\title{Transitive ties}
\usage{
# binary: transitiveties(attr=NULL, levels=NULL)
}
\arguments{
\item{attr}{quantitative attribute (see Specifying Vertex attributes and Levels (\code{?nodal_attributes}) for details.) If set, all three nodes involved ( \eqn{i} , \eqn{j} , and the node on the two-path) must match
on this attribute in order for \eqn{i\rightarrow j}{i-->j} to be counted.}

\item{levels}{TODO (See Specifying Vertex
attributes and Levels (\code{\link[=nodal_attributes]{?nodal_attributes}}) for details.)}
}
\description{
This term adds one statistic, equal to the number of ties
\eqn{i\rightarrow j}{i-->j} such that there exists a two-path from
\eqn{i} to \eqn{j} . (Related to the \code{ttriple} term.)
}
\seealso{
\code{\link{ergmTerm}} for index of model terms currently visible to the package.

\Sexpr[results=rd,stage=render]{ergm:::.formatTermKeywords("ergmTerm", "transitiveties", "subsection")}
}
\concept{categorical nodal attribute}
\concept{directed}
\concept{triad-related}
\concept{undirected}
