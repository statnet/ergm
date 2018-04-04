#  File R/InitWtMHP.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################

InitWtMHP.StdNormal <- function(arguments, nw, response) {
  MHproposal <- list(name = "StdNormal", inputs=NULL)
  MHproposal
}

InitWtMHP.DiscUnif <- function(arguments, nw, response) {
  a <- NVL(arguments$reference$a, -Inf)
  b <- NVL(arguments$reference$b, Inf)
  if(!is.finite(a) || !is.finite(b)) stop('Uniform reference measures that are not bounded are not implemented at this time. Specifiy a and b to be finite.')
  MHproposal <- list(name = "DiscUnif", inputs=c(a,b))
  MHproposal
}

InitWtMHP.DiscUnif2 <- function(arguments, nw, response) {
  a <- NVL(arguments$reference$a, -Inf)
  b <- NVL(arguments$reference$b, Inf)
  if(!is.finite(a) || !is.finite(b)) stop('Uniform reference measures that are not bounded are not implemented at this time. Specifiy a and b to be finite.')
  MHproposal <- list(name = "DiscUnif2", inputs=c(a,b))
  MHproposal
}

InitWtMHP.DiscUnifNonObserved <- function(arguments, nw, response) {
  a <- NVL(arguments$reference$a, -Inf)
  b <- NVL(arguments$reference$b, Inf)
  if(!is.finite(a) || !is.finite(b)) stop('Uniform reference measures that are not bounded are not implemented at this time. Specifiy a and b to be finite.')
  MHproposal <- list(name = "DiscUnifNonObserved", inputs=c(a,b,to_ergm_Cdouble(is.na(nw))))
  MHproposal
}

InitWtMHP.Unif <- function(arguments, nw, response) {
  a <- NVL(arguments$reference$a, -Inf)
  b <- NVL(arguments$reference$b, Inf)
  if(!is.finite(a) || !is.finite(b)) stop('Uniform reference measures that are not bounded are not implemented at this time. Specifiy a and b to be finite.')
  MHproposal <- list(name = "Unif", inputs=c(a,b))
  MHproposal
}

InitWtMHP.UnifNonObserved <- function(arguments, nw, response) {
  a <- NVL(arguments$reference$a, -Inf)
  b <- NVL(arguments$reference$b, Inf)
  if(!is.finite(a) || !is.finite(b)) stop('Uniform reference measures that are not bounded are not implemented at this time. Specifiy a and b to be finite.')
  MHproposal <- list(name = "UnifNonObserved", inputs=c(a,b,to_ergm_Cdouble(is.na(nw))))
  MHproposal
}

InitWtMHP.DistRLE <- function(arguments, nw, response) {
  inputs <- with(arguments$reference,
                 switch(name,
                        Unif = c(0, a, b),
                        DiscUnif = c(1, a, b),
                        StdNormal = c(2, 0, 1),
                        Poisson = c(3, 0),
                        Binomial = c(4, trials, 0.5),
                        Bernoulli = c(4, 1, 0.5)))
  MHproposal <- list(name = "DistRLE", inputs=c(to_ergm_Cdouble(as.rlebdm(arguments$constraints)),inputs), pkgname="ergm")
}
