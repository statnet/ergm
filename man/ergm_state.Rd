% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ergm_state.R, R/summary.ergm_model.R
\docType{data}
\name{ergm_state}
\alias{ergm_state}
\alias{ergm_state.edgelist}
\alias{ergm_state.matrix}
\alias{ergm_state.network}
\alias{is.ergm_state}
\alias{as.edgelist.ergm_state}
\alias{as.matrix.ergm_state}
\alias{as.network.ergm_state_full}
\alias{network.edgecount.ergm_state}
\alias{network.dyadcount.ergm_state_full}
\alias{network.size.ergm_state_full}
\alias{network.naedgecount.ergm_state}
\alias{\%ergmlhs\%.ergm_state_full}
\alias{\%ergmlhs\%<-.ergm_state_full}
\alias{as.rlebdm.ergm_state}
\alias{as.ergm_model.ergm_state_send}
\alias{is.curved.ergm_state_send}
\alias{param_names.ergm_state_send}
\alias{nparam.ergm_state_send}
\alias{update.ergm_state_full}
\alias{ergm_state.ergm_state}
\alias{ERGM_STATE_R_CHANGED}
\alias{ERGM_STATE_C_CHANGED}
\alias{ERGM_STATE_RECONCILED}
\alias{ergm_state_send}
\alias{ergm_state_send.ergm_state_send}
\alias{ergm_state_send.ergm_state_full}
\alias{ergm_state_send.ergm_state_receive}
\alias{update.ergm_state_send}
\alias{ergm_state_receive}
\alias{ergm_state_receive.ergm_state}
\alias{ergm_state_receive.ergm_state_full}
\alias{summary.ergm_state}
\title{A Representation of ERGM state}
\format{
An object of class \code{integer} of length 1.

An object of class \code{integer} of length 1.

An object of class \code{integer} of length 1.
}
\usage{
ergm_state(x, ...)

\method{ergm_state}{edgelist}(
  x,
  nw0,
  model = NULL,
  proposal = NULL,
  stats = NULL,
  ext.state = NULL,
  ...
)

\method{ergm_state}{matrix}(
  x,
  nw0,
  model = NULL,
  proposal = NULL,
  stats = NULL,
  ext.state = NULL,
  ...
)

\method{ergm_state}{network}(x, ...)

is.ergm_state(x)

\method{as.edgelist}{ergm_state}(x, ...)

\method{as.matrix}{ergm_state}(x, matrix.type = NULL, ...)

\method{as.network}{ergm_state_full}(x, ..., populate = TRUE)

\method{network.edgecount}{ergm_state}(x, na.omit = TRUE, ...)

\method{network.dyadcount}{ergm_state_full}(x, na.omit = TRUE, ...)

\method{network.size}{ergm_state_full}(x, ...)

\method{network.naedgecount}{ergm_state}(x, ...)

\method{\%ergmlhs\%}{ergm_state_full}(lhs, setting)

\method{\%ergmlhs\%}{ergm_state_full}(lhs, setting) <- value

\method{as.rlebdm}{ergm_state}(x, ...)

\method{as.ergm_model}{ergm_state_send}(x, ...)

\method{is.curved}{ergm_state_send}(object, ...)

\method{param_names}{ergm_state_send}(object, ...)

\method{nparam}{ergm_state_send}(object, ...)

\method{update}{ergm_state_full}(
  object,
  el = NULL,
  nw0 = NULL,
  model = NULL,
  proposal = NULL,
  stats = NULL,
  ext.state = NULL,
  state = NULL,
  ...
)

\method{ergm_state}{ergm_state}(x, model = NULL, proposal = NULL, stats = NULL, ...)

ERGM_STATE_R_CHANGED

ERGM_STATE_C_CHANGED

ERGM_STATE_RECONCILED

ergm_state_send(x, ...)

\method{ergm_state_send}{ergm_state_send}(x, ...)

\method{ergm_state_send}{ergm_state_full}(x, ...)

\method{ergm_state_send}{ergm_state_receive}(x, ...)

\method{update}{ergm_state_send}(object, state, ...)

ergm_state_receive(x, ...)

\method{ergm_state_receive}{ergm_state}(x, ...)

\method{ergm_state_receive}{ergm_state_full}(x, ...)

\method{summary}{ergm_state}(object, ...)
}
\arguments{
\item{...}{Additional arguments, passed to further methods.}

\item{nw0}{a \code{\link[network:network]{network}} object, whose edges are absent or ignored.}

\item{model}{an \code{\link{ergm_model}} object.}

\item{ext.state}{a list equal to the number of terms in the model,
providing the encoded extended state. This vector is usually
generated by \code{ext.encode()} function of an \code{ergm} term, but it
can be specified directly.}

\item{na.omit}{Whether missing edges should be counted. Note that
missing edge information is not stored.}

\item{state}{An \code{ergm_state} to replace the state with.}
}
\value{
At this time, an \code{ergm_state} object is (subject to change)
a list containing some subset of the following elements, with
\code{el}, \code{ext.state}, and \code{ext.flag} mandatory and others depending
on how it is used:
\describe{

\item{el}{a \code{\link[tibble:tibble]{tibble}} \code{\link[network:as.edgelist]{edgelist}} representing the edge state of the network}

\item{nw0}{a \code{\link[network:network]{network}} object with all edges removed.}

\item{model}{an \code{\link{ergm_model}} object.}

\item{proposal}{an \code{\link{ergm_proposal}} object.}

\item{ext.state}{a list of length equalling to the number of terms in the model.}

\item{ext.flag}{one of \code{ERGM_STATE_R_CHANGED}, \code{ERGM_STATE_C_CHANGED}, and \code{ERGM_STATE_R_RECONCILED}.}

\item{stats}{a numeric vector of network statistics or some other
statistics used to resume.}

\item{uids}{a named list of globally unique ID strings associated with a \code{model} and/or \code{proposal}; for the \code{ergm_state_send} and \code{ergm_state_receive}, these strings may be retained even if these values are set to \code{NULL}}
}
}
\description{
\code{ergm_state} is a family of semi-internal classes for passing
around results of MCMC sampling, particularly when the result is
used to start another MCMC sampler. It is deliberately loosely
specified, and its structure and even name are subject to change.
}
\details{
\code{ergm_state} is actually a hierarchy of classes, defined by what
they can be used for. Specifically,
\describe{

\item{c(\code{ergm_state_receive},\code{ergm_state})}{ needs to contain only \code{el}, \code{ext.state}, and \code{ext.flag}: it is the information that can change in the process of MCMC sampling; it is the one returned by the \verb{*_slave} functions, to minimize the amount of data being sent between nodes in parallel computing.}
\item{c(\code{ergm_state_send},\code{ergm_state_receive},\code{ergm_state})}{ needs the above but also the \code{model} and the \code{proposal}: is needed to initiate MCMC sampling; it is the information required by the \verb{*_slave} functions, again, to minimize the amount of data being sent between nodes in parallel computing.}
\item{c(\code{ergm_state_full}, \code{ergm_state_send},\code{ergm_state_receive},\code{ergm_state})}{ needs the above but also the \code{nw0}: is needed to reconstruct the original network.}

}
}
\section{Methods (by class)}{
\itemize{
\item \code{ergm_state(edgelist)}: a method for constructing an ergm_state from an \code{\link[network:as.edgelist]{edgelist}} object and an empty \code{\link[network:network]{network}}.

\item \code{ergm_state(matrix)}: a method for constructing an ergm_state from a matrix object and an empty \code{\link[network:network]{network}}.

\item \code{ergm_state(network)}: a method for constructing an ergm_state from
a \code{\link[network:network]{network}} object. Note that \code{...} arguments will be passed
directly to the \code{\link[network:as.edgelist]{edgelist}} method.

\item \code{ergm_state(ergm_state)}: a method for constructing an \code{ergm_state}.

}}
\section{Methods (by generic)}{
\itemize{
\item \code{network.edgecount(ergm_state)}: Note that this method fails when
\code{na.omit=FALSE}, since missing edges are not stored.

\item \code{network.naedgecount(ergm_state)}: A stub that returns 0.

\item \code{summary(ergm_state)}: a very low-level function that calculates summary statistics associated with an \code{\link{ergm_state}} object.

}}
\section{Functions}{
\itemize{
\item \code{network.dyadcount(ergm_state_full)}: Note that this method fails with
its default argument, since missing edges are not stored.

\item \code{update(ergm_state_full)}: a method for updating an \code{ergm_state} and reconciling extended state.

}}
\seealso{
\code{\link[=ergm_MCMC_sample]{ergm_MCMC_sample()}} for an example of manually
constructing and manipulating an \code{ergm_state}.
}
\keyword{datasets}
\keyword{internal}
