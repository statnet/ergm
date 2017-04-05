pack.strasnum <- function(s) c(nchar(s), strtoi(charToRaw(s), 16L))
pack.Clistasnum <- function(Clist){
  fnames <- pack.strasnum(Clist$fnamestring)
  snames <- pack.strasnum(Clist$snamestring)
  c(Clist$nterms, fnames, snames, Clist$inputs)
}


## Creates a submodel that does exactly what the model terms passed to
## it would have done.
##
## FIXME: Handle curved terms.

InitErgmTerm.passthrough <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
  f <- a$formula
  if(length(f)==2) f <- nonsimp.update.formula(f, nw~.)
  else nw <- ergm.getnetwork(f)

  m <- ergm.getmodel(f, nw, response=response,...)
  Clist <- ergm.Cprepare(nw, m, response=response)

  inputs <- pack.Clistasnum(Clist)

  gs <- ergm.emptynwstats(m)

  list(name="passthrough_term", coef.names = paste0('passthrough(',m$coef.names,')'), inputs=inputs, dependence=!is.dyad.independent(m), emptynwstats = gs)
}

## Creates a submodel that tracks the given formula.

InitErgmTerm..submodel <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
  f <- a$formula
  if(length(f)==2) f <- nonsimp.update.formula(f, nw~.)
  else nw <- ergm.getnetwork(f)

  m <- ergm.getmodel(f, nw, response=response,...)
  Clist <- ergm.Cprepare(nw, m, response=response)

  inputs <- pack.Clistasnum(Clist)

  gs <- ergm.emptynwstats(m)

  list(name="_submodel_term", coef.names = c(), inputs=inputs, dependence=!is.dyad.independent(m))
}

## Tests .submodel().

InitErgmTerm.submodel.test <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
  f <- a$formula
  if(length(f)==2) f <- nonsimp.update.formula(f, nw~.)
  else nw <- ergm.getnetwork(f)

  m <- ergm.getmodel(f, nw, response=response,...)
  Clist <- ergm.Cprepare(nw, m, response=response)

  gs <- ergm.emptynwstats(m)
  
  list(name="submodel_test_term", coef.names = paste0("submod.test(",m$coef.names,")"), emptynwstats = gs, dependence=!is.dyad.independent(m), auxiliaries = ~.submodel(a$formula))
}

## An auxiliary that exports a double vector that contains the current
## summary statistics of the network for the terms given in the
## formula.

InitErgmTerm..summary <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
  f <- a$formula
  if(length(f)==2) f <- nonsimp.update.formula(f, nw~.)
  else nw <- ergm.getnetwork(f)

  m <- ergm.getmodel(f, nw, response=response,...)
  Clist <- ergm.Cprepare(nw, m, response=response)

  inputs <- pack.Clistasnum(Clist)

  gs <- ergm.emptynwstats(m)

  list(name="_summary_term", coef.names = c(), inputs=c(inputs,gs), dependence=!is.dyad.independent(m))
}


## A term to test (very verbosely) the correctness of the .summary() auxiliary.

InitErgmTerm.summary.test <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL, FALSE),
                      required = c(TRUE, FALSE))

  f <- a$formula
  if(length(f)==2) f <- nonsimp.update.formula(f, nw~.)
  else nw <- ergm.getnetwork(f)

  m <- ergm.getmodel(f, nw, response=response,...)
  
  list(name="summary_test_term", coef.names = 'summ.test', inputs=c(coef.length.model(m)), dependence=!is.dyad.independent(m), auxiliaries=~.summary(a$formula))
}

