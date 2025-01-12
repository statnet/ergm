#  File R/ergm-options.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
#' Global options and term options for the `ergm` package
#'
#' Options set via the built-in [options()] functions that affect `ergm` estimation and options that control the behavior of some terms.
#'
#' @section Global options and defaults:
#'
#' \describe{
#'
#' \item{ergm.eval.loglik = TRUE}{Whether [ergm()] and similar functions will evaluate the likelihood of the fitted model. Can be overridden for a specific call by passing `eval.loglik` argument directly.}
#' 
#' \item{`ergm.loglik.warn_dyads = TRUE`}{Whether log-likelihood evaluation should issue a warning when the effective number of dyads that can vary in the sample space is poorly defined, such as if the degree sequence is constrained.}
#' 
#' \item{`ergm.cluster.retries = 5`}{\CRANpkg{ergm}'s parallel routines implement rudimentary fault-tolerance. This option controls the number of retries for a cluster call before giving up.}
#' 
#' \item{`ergm.term = list()`}{The default term options below.}
#'
#' \item{`ergm.ABI.action = "stop"`}{What to do when \pkg{ergm} detects that one of its extension packages had been compiled with a different version of \pkg{ergm} from the current one that makes changes at the C level that can cause problems. Other choices include \describe{
#'
#' \item{`"stop"`, `"abort"`}{stop with an error}
#'
#' \item{`"warning"`}{warn and proceed}
#'
#' \item{`"message"`, `"inform"`}{print a message and proceed}
#'
#' \item{`"silent"`}{return the value without side-effects}
#'
#' \item{`"disable"`}{skip the check, always returning `TRUE`}
#' }
#' Partial matching is supported. }
#' 
#' }
#' 
#' @section Term options:
#'
#' Term options can be set in three places, in the order of precedence from high to low:
#'
#' 1. As a term argument (not always). For example, `gw.cutoff` below can be set in a \ergmTerm{ergm}{gwesp}{} term by `gwesp(..., cutoff=X)`.
#' 2. For functions such as [`summary`][summary_formula.network] that take `ergm` formulas but do not take a control list, the named arguments passed in as `...`. E.g, `summary(nw~gwesp(.5,fix=TRUE), gw.cutoff=60)` will evaluate the GWESP statistic with its cutoff set to 60.
#' 3. As an element in a `term.options=` list passed via a control function such as [control.ergm()] or, for functions that do not, in a list with that argument name. E.g., `summary(nw~gwesp(.5,fix=TRUE), term.options=list(gw.cutoff=60))` has the same effect.
#' 4. As an element in a global option list `ergm.term` above.
#'
#' The following options are in use by terms in the `ergm` package:
#'
#' \describe{
#'
#' \item{`version`}{A string that can be interpreted as an \R package version. If set, the term will attempt to emulate its behavior as it was that version of `ergm`. Not all past version behaviors are available.}
#' 
#' \item{`gw.cutoff`}{In geometrically weighted terms (`gwesp`, `gwdegree`, etc.) the highest number of shared partners, degrees, etc. for which to compute the statistic. This usually defaults to 30.}
#'
#' \item{`cache.sp`}{Whether the \ergmTerm{ergm}{gwesp}{}, \ergmTerm{ergm}{dgwesp}{}, and similar terms need should use a cache for the dyadwise number of shared partners. This usually improves performance significantly at a modest memory cost, and therefore defaults to `TRUE`, but it can be disabled.}
#'
#' \item{`interact.dependent`}{Whether to allow and how to handle the user attempting to interact dyad-dependent terms (e.g., `absdiff("age"):triangles` or `absdiff("age")*triangles` as opposed to `absdiff("age"):nodefactor("sex")`). Possible values are `"error"` (the default), `"message"`, and `"warning"`, for their respective actions, and `"silent"` for simply processing the term.}
#'
#' }
#'
#' @name ergm-options
#' @aliases term.options ergmTerm-options
NULL
