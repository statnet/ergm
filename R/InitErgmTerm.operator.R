#' Convert a string into a "standard" format for passing as a numeric variable.
#'
#' This function takes a character vector of length 1 and returns a numerical vector encoding it.
#'
#' @param s a character vector of length 1
#'
#' @return a numeric vector concatenating the following:
#' * the number of characters in the input string and
#' * the characters in the input string encoded using ASCII encoding.
#' 
#' This is intended to be decoded by `unpack_str_as_double()` C routines.
#'
#' @export
pack.str_as_num <- function(s) c(nchar(s), strtoi(charToRaw(s), 16L))

#' Serialize model information in a `Clist` into a numeric vector
#'
#' This function takes the output of `ergm.Cprepare` that pertains to
#' the model (not the network) and serializes it by encoding term
#' names and libraries as numeric vectors and concatenating them with
#' inputs and other information.
#'
#' @param Clist output of `ergm.Cprepare`
#'
#' @return a numeric vector concatenating the following:
#' * number of terms in the model;
#' * length of and encoded string of term names;
#' * length of and encoded string of library names; and
#' * vector of intputs to the model.
#' 
#' This is intended to be decoded by `unpack_*Model_as_double()` C routines.
#'
#' @export
pack.Clist_as_num <- function(Clist){
  fnames <- pack.str_as_num(Clist$fnamestring)
  snames <- pack.str_as_num(Clist$snamestring)
  c(Clist$nterms, fnames, snames, Clist$inputs)
}

#' Wrap a submodel's curved specification (if present) for output from an `InitErgmTerm` or `InitWtErgmTerm`.
#'
#' Given a `ergm` model and (optionally) a function with which to wrap
#' parameter names, wrap the calls to its `ergm.eta()` and
#' `ergm.etagrad()` into `map()` and `gradient()` functions, similarly
#' with the `params` element.
#'
#' @param m an `ergm.model` object
#' @param namewrap an optional function taking a character string and
#'   returning a character string, gold on the model's curved
#'   parameter names to wrap them.
#'
#' @return a list with elements `map`, `gradient`, and `params`
#'   suitable for concatenating with an `InitErgmTerm` or
#'   `InitWtErgmTerm` output list.
passthrough.curved.ergm.model <- function(m, namewrap = identity){
  
  if(is.curved(m)){
    map <- function(x, n, ...){
      ergm.eta(x, m$etamap)
    }
    gradient <- function(x, n, ...){
      ergm.etagrad(x, m$etamap)
    }
    params <- rep(list(NULL), coef.length.model(m))
    names(params) <- sapply(coef.names.model(m, canonical=FALSE), namewrap)
  }else map <- gradient <- params <- NULL

  list(map = map, gradient = gradient, params = params)
}

## Creates a submodel that does exactly what the model terms passed to
## it would have done.
##

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

  inputs <- pack.Clist_as_num(Clist)
  
  gs <- ergm.emptynwstats.model(m)
  
  c(list(name="passthrough_term", coef.names = paste0('passthrough(',m$coef.names,')'), inputs=inputs, dependence=!is.dyad.independent(m), emptynwstats = gs),
    passthrough.curved.ergm.model(m, function(x) paste0('passthrough(',x,')')))
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

  inputs <- pack.Clist_as_num(Clist)

  gs <- ergm.emptynwstats.model(m)

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

  gs <- ergm.emptynwstats.model(m)
  
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

  inputs <- pack.Clist_as_num(Clist)

  gs <- ergm.emptynwstats.model(m)

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

