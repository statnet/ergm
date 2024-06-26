% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InitErgmConstraint.hints.R
\name{triadic-ergmHint}
\alias{triadic-ergmHint}
\alias{InitErgmConstraint.triadic}
\alias{triadic-ergmConstraint}
\alias{InitErgmConstraint..triadic}
\alias{.triadic-ergmHint}
\title{Network with strong clustering (triad-closure) effects}
\usage{
# triadic(triFocus = 0.25, type="OTP")

# .triadic(triFocus = 0.25, type = "OTP")
}
\arguments{
\item{triFocus}{A number between 0 and 1, indicating how often triad-focused proposals should be made relative to the standard proposals.}

\item{type}{A string indicating the type of shared partner or path to be considered for directed networks: \code{"OTP"} (default for directed), \code{"ITP"}, \code{"RTP"}, \code{"OSP"}, and \code{"ISP"}; has no effect for undirected. See the section below on Shared partner types for details.}
}
\description{
The network has a high clustering coefficient. This typically results in alternating between the Tie-Non-Tie (TNT) proposal and a triad-focused proposal along the lines of that of \insertCite{WaAt13a;textual}{ergm}.
}
\section{Shared partner types}{
 While there is only one shared partner configuration in the undirected
case, nine distinct configurations are possible for directed graphs, selected
using the \code{type} argument. Currently, terms may be defined with respect to
five of these configurations; they are defined here as follows (using
terminology from Butts (2008) and the \code{relevent} package):
\itemize{
\item Outgoing Two-path (\code{"OTP"}): vertex \eqn{k} is an OTP shared partner of ordered
pair \eqn{(i,j)} iff \eqn{i \to k \to j}{i->k->j}. Also known as "transitive
shared partner".
\item Incoming Two-path (\code{"ITP"}): vertex \eqn{k} is an ITP shared partner of ordered
pair \eqn{(i,j)} iff \eqn{j \to k \to i}{j->k->i}. Also known as "cyclical shared
partner"
\item Reciprocated Two-path (\code{"RTP"}): vertex \eqn{k} is an RTP shared partner of ordered
pair \eqn{(i,j)} iff \eqn{i \leftrightarrow k \leftrightarrow j}{i<->k<->j}.
\item Outgoing Shared Partner (\code{"OSP"}): vertex \eqn{k} is an OSP shared partner of
ordered pair \eqn{(i,j)} iff \eqn{i \to k, j \to k}{i->k, j->k}.
\item Incoming Shared Partner (\code{"ISP"}): vertex \eqn{k} is an ISP shared partner of ordered
pair \eqn{(i,j)} iff \eqn{k \to i, k \to j}{k->i, k->j}.
}

By default, outgoing two-paths (\code{"OTP"}) are calculated. Note that Robins et al. (2009)
define closely related statistics to several of the above, using slightly different terminology.
}

\section{\code{.triadic()} versus \code{triadic()}}{
 If given a bipartite
network, the dotted form will skip silently, whereas the plain
form will raise an error, since triadic effects are not possible
in bipartite networks. The dotted form is thus suitable as a
default argument when the bipartitedness of the network is not
known \emph{a priori}.
}

\references{
\insertAllCited{}
}
\seealso{
\code{\link{ergmHint}} for index of constraints and hints currently visible to the package.

\Sexpr[results=rd,stage=render]{ergm:::.formatTermKeywords("ergmHint", "triadic", "subsection")}
}
\concept{dyad-dependent}
