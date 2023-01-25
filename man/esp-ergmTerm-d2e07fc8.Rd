% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InitErgmTerm.R
\name{esp-ergmTerm}
\alias{esp-ergmTerm}
\alias{InitErgmTerm.esp}
\title{Edgewise shared partners}
\usage{
# binary: esp(d)
}
\arguments{
\item{d}{a vector of distinct integers}
}
\description{
This is just like the \code{dsp} term, except this term adds one network
statistic to the model for each element in \code{d} where the \eqn{i} th such
statistic equals the number of edges (rather than dyads) in the
network with exactly \code{d[i]} shared partners. This term can be used with
directed and undirected networks.
}
\note{
This term takes an additional term option (see
\code{\link[=ergm-options]{options?ergm}}), \code{cache.sp}, controlling whether
the implementation will cache the number of shared partners for
each dyad in the network; this is usually enabled by default.

For directed networks, only outgoing two-path ("OTP") shared partners
are counted.  In other words, for a (directed) edge \code{i -> j} in a directed graph, the number of
shared partners counted by \code{esp} is the number of nodes \code{k} that have edges
\code{i -> k -> j}.  (These may also be called homogeneous shared partners.)  To
count other types of shared partners instead, see \code{desp}.
}
\seealso{
\code{\link{ergmTerm}} for index of model terms currently visible to the package.

\Sexpr[results=rd,stage=render]{ergm:::.formatTermKeywords("ergmTerm", "esp", "subsection")}
}
\concept{directed}
\concept{undirected}