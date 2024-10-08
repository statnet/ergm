% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ergm_model.R, R/is.curved.R,
%   R/is.dyad.independent.R, R/nparam.R, R/param_names.R
\name{ergm_model}
\alias{ergm_model}
\alias{c.ergm_model}
\alias{as.ergm_model}
\alias{as.ergm_model.ergm_model}
\alias{as.ergm_model.formula}
\alias{is.curved.ergm_model}
\alias{is.dyad.independent.ergm_model}
\alias{nparam.ergm_model}
\alias{param_names.ergm_model}
\alias{param_names<-.ergm_model}
\title{Internal representation of an \code{ergm} network model}
\usage{
ergm_model(
  formula,
  nw = NULL,
  silent = FALSE,
  ...,
  term.options = list(),
  extra.aux = list(),
  env = globalenv(),
  offset.decorate = TRUE,
  terms.only = FALSE
)

\method{c}{ergm_model}(...)

as.ergm_model(x, ...)

\method{as.ergm_model}{ergm_model}(x, ...)

\method{as.ergm_model}{formula}(x, ...)

\method{is.curved}{ergm_model}(object, ...)

\method{is.dyad.independent}{ergm_model}(object, ..., ignore_aux = TRUE)

\method{nparam}{ergm_model}(object, canonical = FALSE, offset = NA, byterm = FALSE, ...)

\method{param_names}{ergm_model}(object, canonical = FALSE, offset = NA, ...)

\method{param_names}{ergm_model}(object, canonical = FALSE, ...) <- value
}
\arguments{
\item{formula}{An \code{\link[=ergm]{ergm()}}
formula of the form \code{network ~ model.term(s)} or \code{~
model.term(s)} or a \code{\link[statnet.common:term_list]{term_list}} object, typically constructed from a formula's LHS.}

\item{nw}{The network of interest, optionally instrumented with \code{\link[=ergm_preprocess_response]{ergm_preprocess_response()}} to have a response attribute specification; if passed, the LHS of \code{formula} is ignored. This is the recommended usage.}

\item{silent}{logical, whether to print the warning messages from the
initialization of each model term.}

\item{\dots}{additional parameters for model formulation}

\item{term.options}{A list of additional arguments to be passed to term initializers. See \code{\link[=term.options]{? term.options}}.}

\item{extra.aux}{a list of auxiliary request formulas required elsewhere; if named, the resulting \code{slots.extra.aux} will also be named.}

\item{env}{a throwaway argument needed to prevent conflicts with some usages of \code{ergm_model}. The initialization environment is \emph{always} taken from the \code{formula}.}

\item{offset.decorate}{logical; whether offset coefficient and parameter names should be enclosed in \code{"offset()"}.}

\item{terms.only}{logical; whether auxiliaries, eta map, and UID constructions should be skipped. This is useful for submodels.}

\item{x}{object to be converted to an \code{ergm_model}.}

\item{object}{An \code{ergm_model} object.}

\item{ignore_aux}{A flag to specify whether a dyad-dependent
auxiliary should make the model dyad-dependent or should be
ignored.}

\item{canonical}{Whether the canonical (eta) parameters or the curved (theta) parameters are used.}

\item{offset}{If \code{NA} (the default), all model terms are counted;
if \code{TRUE}, only offset terms are counted; and if
\code{FALSE}, offset terms are skipped.}

\item{byterm}{Whether to return a vector of the numbers of
coefficients for each term.}

\item{value}{For \code{\link[=param_names<-]{param_names<-()}}, either a character vector equal
in length to the number of parameters of the specified type
(though recycled as needed), or a \code{\link{list}} of two character
vectors, one for non-canonical, the other for canonical, in which
case \verb{canonical=} will be ignored. \code{NA} elements preserve
existing name.}
}
\value{
\code{ergm_model} returns an  \code{ergm_model} object as a list
containing:
\item{terms}{a list of terms and 'term components' initialized by the
appropriate \code{InitErgmTerm.X} function.}
\item{etamap}{the theta -> eta mapping as a list returned from
<ergm.etamap>}
\item{uid}{a string generated with the model, \UIDalgo; different models are, generally, guaranteed to have different strings, but identical models are not guaranteed to have the same string}
}
\description{
These methods are generally not called directly by users, but may
be employed by other depending packages.
\code{ergm_model} constructs it from a formula or a term list. Each term is
initialized via the \code{InitErgmTerm} functions to create a
\code{ergm_model} object.
}
\section{Methods (by generic)}{
\itemize{
\item \code{c(ergm_model)}: A method for concatenating terms of two or more initialized models.

\item \code{is.curved(ergm_model)}: Tests whether the model is curved.

\item \code{is.dyad.independent(ergm_model)}: Tests whether the model is dyad-independent.

\item \code{nparam(ergm_model)}: Number of parameters of the model.

\item \code{param_names(ergm_model)}: Parameter names of the model.

\item \code{param_names(ergm_model) <- value}: Rename the parameters.

}}
\note{
This API is not to be considered fixed and may change between versions. However, an effort will be made to ensure that the methods of this class remain stable.

Earlier versions also had an optional \verb{response=} parameter that, if not \code{NULL}, switched to valued mode and used the edge attribute named in \verb{response=} as the response. This is no longer used; instead, the response is to be set on \code{nw} via \code{ergm_preprocess_response(nw, response)}.
}
\seealso{
\code{\link[=summary.ergm_model]{summary.ergm_model()}}, \code{\link[=ergm_preprocess_response]{ergm_preprocess_response()}}
}
\keyword{internal}
