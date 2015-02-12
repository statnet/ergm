#  File R/InitWtMHP.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
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

InitWtMHP.DiscUnifNonObserved <- function(arguments, nw, response) {
  a <- NVL(arguments$reference$a, -Inf)
  b <- NVL(arguments$reference$b, Inf)
  if(!is.finite(a) || !is.finite(b)) stop('Uniform reference measures that are not bounded are not implemented at this time. Specifiy a and b to be finite.')
  MHproposal <- list(name = "DiscUnifNonObserved", inputs=c(a,b,ergm.Cprepare.miss(nw)))
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
  MHproposal <- list(name = "UnifNonObserved", inputs=c(a,b,ergm.Cprepare.miss(nw)))
  MHproposal
}
