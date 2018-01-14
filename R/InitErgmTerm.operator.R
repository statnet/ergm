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
#' @param m an `ergm_model` object
#' @param namewrap an optional function taking a character string and
#'   returning a character string, gold on the model's curved
#'   parameter names to wrap them.
#'
#' @return a list with elements `map`, `gradient`, and `params`
#'   suitable for concatenating with an `InitErgmTerm` or
#'   `InitWtErgmTerm` output list.
passthrough.curved.ergm_model <- function(m, namewrap = identity){
  
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
    passthrough.curved.ergm_model(m, function(x) paste0('passthrough(',x,')')))
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

InitErgmTerm.F <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "form"),
                      vartypes = c("formula", "formula"),
                      defaultvalues = list(NULL, "sum"),
                      required = c(TRUE, FALSE))
  form <- a$form
  f <- a$formula
  if(length(f)==2) f <- nonsimp.update.formula(f, nw~.)
  else nw <- ergm.getnetwork(f)

  m <- ergm.getmodel(f, nw,...)
  Clist <- ergm.Cprepare(nw, m)
  inputs <- pack.Clist_as_num(Clist)
  
  gs <- ergm.emptynwstats.model(m)

  form.name <- deparse(form[[length(form)]])
  name <- "filter_term_form"
  auxiliaries <- ~.filter.formula.net(form)
  
  c(list(name=name,
         coef.names = paste0(form.name,'(',m$coef.names,')'),
         inputs=inputs,
         dependence=!is.dyad.independent(m),
         emptynwstats = gs,
         auxiliaries=auxiliaries),
    passthrough.curved.ergm_model(m, function(x) paste0(form.name,'(',x,')')))
}

InitErgmTerm..filter.formula.net <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  # Form is a model.
  f<-a$formula
  if(length(f)==2) f <- nonsimp.update.formula(f, nw~.)
  else nw <- ergm.getnetwork(f)
  
  m <- ergm.getmodel(f, nw, response=response,...)

  if(!is.dyad.independent(m) || coef.length.model(m)!=1) stop("The filter test formula must be dyad-independent and have exactly one statistc.")

  Clist <- ergm.Cprepare(nw, m)
  inputs <- pack.Clist_as_num(Clist)

  gs <- ergm.emptynwstats.model(m)
  if(gs!=0) stop("At this time, the filter test term must have the property that its dyadwise components are 0 for 0-valued relations. This limitation may be removed in the future.")
  
  list(name="_filter_formula_net", inputs=c(inputs), depenence=FALSE)
}

InitErgmTerm.Offset <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "coef"),
                      vartypes = c("formula", "numeric"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, TRUE))
  f <- a$formula
  if(length(f)==2) f <- nonsimp.update.formula(f, nw~.)
  else nw <- ergm.getnetwork(f)

  m <- ergm.getmodel(f, nw, response=response,...)
  coef <- rep(a$coef, length(m$coef.names))
  Clist <- ergm.Cprepare(nw, m, response=response)

  inputs <- pack.Clist_as_num(Clist)
  
  gs <- ergm.emptynwstats.model(m)
  
  list(name="passthrough_term", coef.names = paste0('Offset(',m$coef.names,',',coef,')'), inputs=inputs, dependence=!is.dyad.independent(m), emptynwstats = gs,
       params=list(),
       map = function(x, n, ...){
         ergm.eta(coef, m$etamap)
       },
       gradient = function(x, n, ...){
         matrix(NA, 0, length(coef))
       }
       )
}

#' Return a symmetrized version of a binary network
#'
#' @param x an object representing a network.
#' @param rule a string specifying how the network is to be
#'   symmetrized; see [sna::symmetrize()] for details.
#' @param ... additional arguments to [sna::symmetrize()].
#' @export
symmetrize <- function(x, rule=c("weak","strong","upper","lower"), ...){
  UseMethod("symmetrize")
}

#' @describeIn symmetrize
#'
#' The default method, passing the input on to [sna::symmetrize()].
#' 
#' @export
symmetrize.default <- function(x, rule=c("weak","strong","upper","lower"), ...){
  sna::symmetrize(x, rule=rule, ...)
}

#' @describeIn symmetrize
#'
#' A method for [`network`] objects, which preserves network and vertex attributes.
#' 
#' @export
symmetrize.network <- function(x, rule=c("weak","strong","upper","lower"), ...){
  rule <- match.arg(rule)
  el <- sna::symmetrize(x, rule=rule, return.as.edgelist=TRUE, ...)
  o <- network.initialize(network.size(x), directed=FALSE, bipartite=x%n%"bipartite", loops=has.loops(x), hyper=is.hyper(x), multiple=is.multiplex(x))
  el <- el[seq_len(nrow(el))/2,-3,drop=FALSE]
  o <- network.edgelist(el, o)
  nvattr.copy.network(o, x)
}

InitErgmTerm.Undir <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist, directed = TRUE,
                      varnames = c("formula", "rule"),
                      vartypes = c("formula", "character"),
                      defaultvalues = list(NULL, "weak"),
                      required = c(TRUE, FALSE))
  RULES <- c("weak","strong","upper","lower")
  rule <- match.arg(a$rule, RULES)

  f <- a$formula
  if(length(f)==3) nw <- ergm.getnetwork(f)
  if(is.directed(nw)) nw <- symmetrize(nw, rule)
  
  if(length(f)==2) f <- nonsimp.update.formula(f, nw~.)

  m <- ergm.getmodel(f, nw,...)
  Clist <- ergm.Cprepare(nw, m)
  inputs <- pack.Clist_as_num(Clist)
  
  gs <- ergm.emptynwstats.model(m)

  auxiliaries <- ~.undir.net(rule)
  
  c(list(name="undir",
         coef.names = paste0('Undir(',m$coef.names,')'),
         inputs=c(which(RULES==rule),inputs),
         dependence=!is.dyad.independent(m) || rule%in%c("weak","strong"),
         emptynwstats = gs,
         auxiliaries=auxiliaries),
    passthrough.curved.ergm_model(m, function(x) paste0('Undir(',x,')')))
}

