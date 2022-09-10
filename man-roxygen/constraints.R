#  File man-roxygen/constraints.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
#' @param constraints A formula specifying one or more constraints
#' on the support of the distribution of the networks being modeled. Multiple constraints
#' may be given, separated by \dQuote{+} and \dQuote{-} operators. See
#' [`ergmConstraint`] for the detailed explanation of
#' their semantics and also for an indexed list of the constraints visible to the \CRANpkg{ergm} package.
#'
#' The default is to have no constraints except those provided through
#' the [`ergmlhs`] API.
#'
#' Together with the model terms in the formula and the reference measure, the constraints
#' define the distribution of networks being modeled.
#' 
#' It is also possible to specify a proposal function directly either
#' by passing a string with the function's name (in which case,
#' arguments to the proposal should be specified through the
#' \code{MCMC.prop.args} argument to the relevant control function, or
#' by giving it on the LHS of the hints formula to `MCMC.prop`
#' argument to the control function. This will override
#' the one chosen automatically.
#' 
#' Note that not all possible combinations of constraints and reference
#' measures are supported. However, for relatively simple constraints
#' (i.e., those that simply permit or forbid specific dyads or sets of
#' dyads from changing), arbitrary combinations should be possible.
#' 
