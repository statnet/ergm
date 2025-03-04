% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InitErgmTerm.operator.R,
%   R/InitWtErgmTerm.operator.R
\name{Symmetrize-ergmTerm}
\alias{Symmetrize-ergmTerm}
\alias{InitErgmTerm.Symmetrize}
\alias{InitWtErgmTerm.Symmetrize}
\title{Evaluation on symmetrized (undirected) network}
\usage{
# binary: Symmetrize(formula, rule="weak")

# valued: Symmetrize(formula, rule="weak")
}
\arguments{
\item{formula}{a one-sided \code{\link[=ergm]{ergm()}}-style formula with the terms to be evaluated}

\item{rule}{one of \code{"weak"}, \code{"strong"}, \code{"upper"}, \code{"lower"}, \code{"max"}, or \code{"min"}}
}
\description{
Evaluates the terms in \code{formula} on an undirected network
constructed by symmetrizing the LHS network using one of four rules:
\enumerate{
\item \code{"weak"}/\code{"max"} A tie \eqn{(i,j)} is present in the constructed
network if the LHS network has either tie \eqn{(i,j)} or
\eqn{(j,i)} (or both). For a valued network, the maximum is used.
\item \code{"strong"}/\code{"min"} A tie \eqn{(i,j)} is present in the constructed
network if the LHS network has both tie \eqn{(i,j)} and tie
\eqn{(j,i)}. For a valued network, the minimum is used.
\item "upper" A tie \eqn{(i,j)} is present in the constructed
network if the LHS network has tie \eqn{(\min(i,j),\max(i,j))}{(min(i,j), max(i,j))} :
the upper triangle of the LHS network.
\item "lower" A tie \eqn{(i,j)} is present in the constructed
network if the LHS network has tie \eqn{(\max(i,j),\min(i,j))}{(max(i,j), min(i,j))} :
the lower triangle of the LHS network.
}
}
\seealso{
\code{\link{ergmTerm}} for index of model terms currently visible to the package.

\Sexpr[results=rd,stage=render]{ergm:::.formatTermKeywords("ergmTerm", "Symmetrize", "subsection")}
}
\concept{directed}
\concept{operator}
