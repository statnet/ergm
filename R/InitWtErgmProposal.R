#  File R/InitWtErgmProposal.R in package ergm, part of the Statnet suite of
#  packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################

#' @templateVar name DiscUnif2
#' @aliases InitWtErgmProposal.DiscUnif2
#' @title Discrete Uniform proposal that proposes changes to two dyads
#' @description Mainly used for testing
#' @template ergmProposal-general
NULL
InitWtErgmProposal.DiscUnif2 <- function(arguments, nw) {
  a <- NVL(arguments$reference$arguments$a, -Inf)
  b <- NVL(arguments$reference$arguments$b, Inf)
  if(!is.finite(a) || !is.finite(b)) stop('Uniform reference measures that are not bounded are not implemented at this time. Specifiy a and b to be finite.')
  proposal <- list(name = "DiscUnif2", inputs=c(a,b))
  proposal
}


#' @templateVar name Dist
#' @aliases InitWtErgmProposal.Dist
#' @title Random toggle for a number of reference distributions
#' @description Implements the `DyadGen` API and the
#'   \ergmConstraint{ergm}{ChangeStats}{()} API.
#' @template ergmProposal-general
NULL
InitWtErgmProposal.Dist <- function(arguments, nw) {
  iinputs <- c(Unif = 0,
               DiscUnif = 1,
               StdNormal = 2,
               Poisson = 3,
               Binomial = 4,
               Bernoulli = 4)[arguments$reference$name]
  inputs <- with(arguments$reference$arguments,
                 switch(arguments$reference$name,
                        Unif = c(a, b),
                        DiscUnif = c(a, b),
                        StdNormal = c(NVL(arguments$sd, 0.2)),
                        Binomial = c(trials),
                        Bernoulli = c(1)))
  proposal <- c(list(name = "Dist", dyadgen = ergm_dyadgen_select(arguments, nw),
                     iinputs = iinputs, inputs = inputs, pkgname = "ergm"),
                ergm_constrain_changestats(arguments))
}
