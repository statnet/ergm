#  File R/InitErgmTerm.projection.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
#' @templateVar name Project
#' @title Evaluation on a projection of a bipartite network
#'
#' @description This operator on a bipartite network evaluates the
#'   formula on the undirected, valued network constructed by
#'   projecting it onto its specified mode. `Proj1(formula)` and
#'   `Proj2(formula)` are aliases for `Project(formula, 1)` and
#'   `Project(formula, 2)`, respectively.
#'
#' @usage
#' # binary: Project(formula, mode)
#' @template ergmTerm-formula
#' @param mode the mode onto which to project: 1 or 2
#'
#' @template ergmTerm-general
#'
#' @concept operator
#' @concept bipartite
InitErgmTerm.Project <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist, bipartite = TRUE,
                      varnames = c("formula", "mode"),
                      vartypes = c("formula", "numeric"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, TRUE))

  bip <- as.integer(nw %n% "bipartite")
  n <- as.integer(network.size(nw))

  mode <- a$mode
  if(! mode %in% 1:2) ergm_Init_stop(sQuote("mode"), " must be 1 or 2.")

  ### Construct an empty network with the correct structure.
  pnw <- nw
  pnw %n% "bipartite" <- FALSE
  pnw[,] <- 0
  pnw <- switch(mode,
                get.inducedSubgraph(pnw, seq_len(bip)),
                get.inducedSubgraph(pnw, bip + seq_len(n-bip)))
  pnw %ergmlhs% "response" <- structure("w", valued = TRUE)

  m <- ergm_model(a$formula, pnw, ..., offset.decorate=FALSE)
  ergm_no_ext.encode(m)

  auxiliaries <- trim_env(~.project.net(mode), "mode")

  c(list(name="on_proj_net", iinputs = mode,
         submodel = m, dependence = TRUE,
         auxiliaries=auxiliaries),
    wrap.ergm_model(m, pnw, ergm_mk_std_op_namewrap(paste0('Proj', mode))))
}

#' @templateVar name Project
#' @template ergmTerm-rdname
#' @aliases Proj1-ergmTerm
#' @usage
#' # binary: Proj1(formula)
InitErgmTerm.Proj1 <- function(nw, arglist, ...){
  arglist[["mode"]] <- 1L
  f <- InitErgmTerm.Project
  f(nw, arglist, ...)
}

#' @templateVar name Project
#' @template ergmTerm-rdname
#' @aliases Proj2-ergmTerm
#' @usage
#' # binary: Proj2(formula)
InitErgmTerm.Proj2 <- function(nw, arglist, ...){
  arglist[["mode"]] <- 2L
  f <- InitErgmTerm.Project
  f(nw, arglist, ...)
}

InitErgmTerm..project.net <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist, bipartite = TRUE,
                      varnames = c("mode"),
                      vartypes = c("numeric"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  mode <- a$mode
  if(! mode %in% 1:2) ergm_Init_stop(sQuote("mode"), " must be 1 or 2.")
  list(name=paste0("_proj_net"), iinputs = mode, coef.names = c(), dependence=TRUE)
}
