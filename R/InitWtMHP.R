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
