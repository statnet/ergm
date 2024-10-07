#  File R/InitErgmReference.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################

#' @templateVar name Bernoulli
#' @title Bernoulli reference
#' @description Specifies each
#'   dyad's baseline distribution to be Bernoulli with probability of
#'   the tie being \eqn{0.5} . This is the only reference measure used
#'   in binary mode.
#'
#' @usage
#' # Bernoulli
#'
#' @template ergmReference-general
#'
#' @concept discrete
#' @concept finite
#' @concept binary
#' @concept nonnegative
InitErgmReference.Bernoulli <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist)
  list(name="Bernoulli", init_methods=c("MPLE", "CD", "zeros"))
}

#' @templateVar name StdNormal
#' @title Standard Normal reference
#' @description Specifies each dyad's baseline distribution to be the normal distribution
#'   with mean 0 and variance 1.
#'
#' @usage
#' # StdNormal
#'
#' @template ergmReference-general
#'
#' @concept continuous
InitErgmReference.StdNormal <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist)
  list(name="StdNormal", init_methods=c("CD","zeros"))
}

#' @templateVar name Unif
#' @title Continuous Uniform reference
#' @description Specifies each dyad's baseline distribution to be continuous uniform
#'   between `a` and `b`: \eqn{h(y)=1} , with the support being `[a, b]`.
#'
#' @usage
#' # Unif(a,b)
#' @param a,b minimum and maximum to the baseline discrete uniform distribution, both inclusive. Both values must be finite.
#'
#' @template ergmReference-general
#'
#' @concept continuous
InitErgmReference.Unif <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("a", "b"),
                      vartypes = c("numeric", "numeric"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, TRUE))
  list(name="Unif", arguments=list(a=a$a, b=a$b), init_methods=c("CD","zeros"))
}

#' @templateVar name DiscUnif
#' @title Discrete Uniform reference
#' @description Specifies each dyad's baseline distribution to be discrete uniform
#'   between `a` and `b` (both inclusive): \eqn{h(y)=1} , with
#'   the support being
#'   `a`, `a+1`, \ldots, `b-1`, `b`.
#'
#' @usage
#' # DiscUnif(a,b)
#' @param a,b minimum and maximum to the baseline discrete uniform distribution, both inclusive. Both values must be finite.
#'
#' @template ergmReference-general
#'
#' @concept discrete
#' @concept finite
InitErgmReference.DiscUnif <- function(nw, arglist, a, b, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("a", "b"),
                      vartypes = c("numeric", "numeric"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, TRUE))
  if(a$a!=round(a$a) || a$b != round(a$b)) ergm_Init_stop("arguments ", sQuote("a"), " and ", sQuote("b"), "must be integers")
  list(name="DiscUnif", arguments=list(a=a$a, b=a$b), init_methods=c("CD","zeros"))
}
