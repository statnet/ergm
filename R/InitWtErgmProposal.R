#  File R/InitWtErgmProposal.R in package ergm, part of the Statnet suite of
#  packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################

#' @templateVar name StdNormal
#' @aliases InitWtErgmProposal.StdNormal
#' @title TODO
#' @description TODO
#' @template ergmProposal-general
NULL
InitWtErgmProposal.StdNormal <- function(arguments, nw) {
  proposal <- list(name = "StdNormal", inputs=NULL)
  proposal
}

#' @templateVar name DiscUnif
#' @aliases InitWtErgmProposal.DiscUnif
#' @title TODO
#' @description TODO
#' @template ergmProposal-general
NULL
InitWtErgmProposal.DiscUnif <- function(arguments, nw) {
  a <- NVL(arguments$reference$arguments$a, -Inf)
  b <- NVL(arguments$reference$arguments$b, Inf)
  if(!is.finite(a) || !is.finite(b)) ergm_Init_stop('Uniform reference measures that are not bounded are not implemented at this time. Specifiy a and b to be finite.')
  proposal <- list(name = "DiscUnif", inputs=c(a,b))
  proposal
}

#' @templateVar name DiscUnif2
#' @aliases InitWtErgmProposal.DiscUnif2
#' @title TODO
#' @description TODO
#' @template ergmProposal-general
NULL
InitWtErgmProposal.DiscUnif2 <- function(arguments, nw) {
  a <- NVL(arguments$reference$arguments$a, -Inf)
  b <- NVL(arguments$reference$arguments$b, Inf)
  if(!is.finite(a) || !is.finite(b)) stop('Uniform reference measures that are not bounded are not implemented at this time. Specifiy a and b to be finite.')
  proposal <- list(name = "DiscUnif2", inputs=c(a,b))
  proposal
}

#' @templateVar name DiscUnifNonObserved
#' @aliases InitWtErgmProposal.DiscUnifNonObserved
#' @title TODO
#' @description TODO
#' @template ergmProposal-general
NULL
InitWtErgmProposal.DiscUnifNonObserved <- function(arguments, nw) {
  a <- NVL(arguments$reference$arguments$a, -Inf)
  b <- NVL(arguments$reference$arguments$b, Inf)
  if(!is.finite(a) || !is.finite(b)) ergm_Init_stop('Uniform reference measures that are not bounded are not implemented at this time. Specifiy a and b to be finite.')
  proposal <- list(name = "DiscUnifNonObserved", inputs=c(a,b,to_ergm_Cdouble(is.na(nw))))
  proposal
}

#' @templateVar name Unif
#' @aliases InitWtErgmProposal.Unif
#' @title TODO
#' @description TODO
#' @template ergmProposal-general
NULL
InitWtErgmProposal.Unif <- function(arguments, nw) {
  a <- NVL(arguments$reference$arguments$a, -Inf)
  b <- NVL(arguments$reference$arguments$b, Inf)
  if(!is.finite(a) || !is.finite(b)) ergm_Init_stop('Uniform reference measures that are not bounded are not implemented at this time. Specifiy a and b to be finite.')
  proposal <- list(name = "Unif", inputs=c(a,b))
  proposal
}

#' @templateVar name UnifNonObserved
#' @aliases InitWtErgmProposal.UnifNonObserved
#' @title TODO
#' @description TODO
#' @template ergmProposal-general
NULL
InitWtErgmProposal.UnifNonObserved <- function(arguments, nw) {
  a <- NVL(arguments$reference$arguments$a, -Inf)
  b <- NVL(arguments$reference$arguments$b, Inf)
  if(!is.finite(a) || !is.finite(b)) ergm_Init_stop('Uniform reference measures that are not bounded are not implemented at this time. Specifiy a and b to be finite.')
  proposal <- list(name = "UnifNonObserved", inputs=c(a,b,to_ergm_Cdouble(is.na(nw))))
  proposal
}

#' @templateVar name DistRLE
#' @aliases InitWtErgmProposal.DistRLE
#' @title TODO
#' @description TODO
#' @template ergmProposal-general
NULL
InitWtErgmProposal.DistRLE <- function(arguments, nw) {
  inputs <- with(arguments$reference$arguments,
                 switch(arguments$reference$name,
                        Unif = c(0, a, b),
                        DiscUnif = c(1, a, b),
                        StdNormal = c(2, 0, 1),
                        Poisson = c(3, 1),
                        Binomial = c(4, trials, 0.5),
                        Bernoulli = c(4, 1, 0.5)))
  proposal <- list(name = "DistRLE", inputs=c(to_ergm_Cdouble(as.rlebdm(arguments$constraints)),inputs), pkgname="ergm")
}
