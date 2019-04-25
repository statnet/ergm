#  File R/InitWtErgmProposal.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################

InitWtErgmProposal.StdNormal <- function(arguments, nw, response) {
  proposal <- list(name = "StdNormal", inputs=NULL)
  proposal
}

InitWtErgmProposal.DiscUnif <- function(arguments, nw, response) {
  a <- NVL(arguments$reference$arguments$a, -Inf)
  b <- NVL(arguments$reference$arguments$b, Inf)
  if(!is.finite(a) || !is.finite(b)) ergm_Init_abort('Uniform reference measures that are not bounded are not implemented at this time. Specifiy a and b to be finite.')
  proposal <- list(name = "DiscUnif", inputs=c(a,b))
  proposal
}

InitWtErgmProposal.DiscUnifNonObserved <- function(arguments, nw, response) {
  a <- NVL(arguments$reference$arguments$a, -Inf)
  b <- NVL(arguments$reference$arguments$b, Inf)
  if(!is.finite(a) || !is.finite(b)) ergm_Init_abort('Uniform reference measures that are not bounded are not implemented at this time. Specifiy a and b to be finite.')
  proposal <- list(name = "DiscUnifNonObserved", inputs=c(a,b,to_ergm_Cdouble(is.na(nw))))
  proposal
}

InitWtErgmProposal.Unif <- function(arguments, nw, response) {
  a <- NVL(arguments$reference$arguments$a, -Inf)
  b <- NVL(arguments$reference$arguments$b, Inf)
  if(!is.finite(a) || !is.finite(b)) ergm_Init_abort('Uniform reference measures that are not bounded are not implemented at this time. Specifiy a and b to be finite.')
  proposal <- list(name = "Unif", inputs=c(a,b))
  proposal
}

InitWtErgmProposal.UnifNonObserved <- function(arguments, nw, response) {
  a <- NVL(arguments$reference$arguments$a, -Inf)
  b <- NVL(arguments$reference$arguments$b, Inf)
  if(!is.finite(a) || !is.finite(b)) ergm_Init_abort('Uniform reference measures that are not bounded are not implemented at this time. Specifiy a and b to be finite.')
  proposal <- list(name = "UnifNonObserved", inputs=c(a,b,to_ergm_Cdouble(is.na(nw))))
  proposal
}

InitWtErgmProposal.DistRLE <- function(arguments, nw, response) {
  inputs <- with(arguments$reference$arguments,
                 switch(arguments$reference$name,
                        Unif = c(0, a, b),
                        DiscUnif = c(1, a, b),
                        StdNormal = c(2, 0, 1),
                        Poisson = c(3, 0),
                        Binomial = c(4, trials, 0.5),
                        Bernoulli = c(4, 1, 0.5)))
  proposal <- list(name = "DistRLE", inputs=c(to_ergm_Cdouble(as.rlebdm(arguments$constraints)),inputs), pkgname="ergm")
}
