#  File R/InitErgmTerm.blockop.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2023 Statnet Commons
################################################################################
## Creates a submodel that ignores any edges not within the
## blocks.

#' @templateVar name NodematchFilter
#' @title Filtering on nodematch
#' @description Evaluates the terms specified in `formula` on a network
#'   constructed by taking \eqn{y} and removing any edges for which
#'   `attrname(i)!=attrname(j)` .
#'
#' @usage
#' # binary: NodematchFilter(formula, attrname)
#' @param formula formula to be evaluated
#' @param attrname a character vector giving one or more names of attributes in the network's vertex attribute list.
#'
#' @template ergmTerm-general
#'
#' @concept operator
InitErgmTerm.NodematchFilter <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "attrname"),
                      vartypes = c("formula", "character"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, TRUE))

  m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)

  attrname <- a$attrname
  c(list(name="on_blockdiag_net", submodel=m, auxiliaries=trim_env(~.blockdiag.net(attrname), "attrname")),
    wrap.ergm_model(m, nw, ergm_mk_std_op_namewrap('NodematchFilter',a$attrname)))
}

