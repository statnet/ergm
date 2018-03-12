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
#' @export passthrough.curved.ergm_model
passthrough.curved.ergm_model <- function(m, namewrap = identity){
  
  if(is.curved(m)){
    map <- function(x, n, ...){
      ergm.eta(x, m$etamap)
    }
    gradient <- function(x, n, ...){
      ergm.etagrad(x, m$etamap)
    }
    params <- rep(list(NULL), nparam(m))
    names(params) <- sapply(param_names(m, canonical=FALSE), namewrap)
  }else map <- gradient <- params <- NULL

  list(map = map, gradient = gradient, params = params, minpar=m$etamap$mintheta, maxpar=m$etamap$maxtheta)
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
  if(length(f)==2) f <- nonsimp_update.formula(f, nw~.)
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
  if(length(f)==2) f <- nonsimp_update.formula(f, nw~.)
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
  if(length(f)==2) f <- nonsimp_update.formula(f, nw~.)
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
  if(length(f)==2) f <- nonsimp_update.formula(f, nw~.)
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
  if(length(f)==2) f <- nonsimp_update.formula(f, nw~.)
  else nw <- ergm.getnetwork(f)

  m <- ergm.getmodel(f, nw, response=response,...)
  
  list(name="summary_test_term", coef.names = 'summ.test', inputs=c(nparam(m)), dependence=!is.dyad.independent(m), auxiliaries=~.summary(a$formula))
}

InitErgmTerm.F <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "form"),
                      vartypes = c("formula", "formula"),
                      defaultvalues = list(NULL, "sum"),
                      required = c(TRUE, FALSE))
  form <- a$form
  f <- a$formula
  if(length(f)==2) f <- nonsimp_update.formula(f, nw~.)
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
  if(length(f)==2) f <- nonsimp_update.formula(f, nw~.)
  else nw <- ergm.getnetwork(f)
  
  m <- ergm.getmodel(f, nw, response=response,...)

  if(!is.dyad.independent(m) || nparam(m)!=1) stop("The filter test formula must be dyad-independent and have exactly one statistc.")

  Clist <- ergm.Cprepare(nw, m)
  inputs <- pack.Clist_as_num(Clist)

  gs <- ergm.emptynwstats.model(m)
  if(gs!=0) stop("At this time, the filter test term must have the property that its dyadwise components are 0 for 0-valued relations. This limitation may be removed in the future.")
  
  list(name="_filter_formula_net", inputs=c(inputs), depenence=FALSE)
}

InitErgmTerm.Offset <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "coef", "which"),
                      vartypes = c("formula", "numeric", "logical,numeric,character"),
                      defaultvalues = list(NULL, 0, TRUE),
                      required = c(TRUE, FALSE, FALSE))
  f <- a$formula
  if(length(f)==2) f <- nonsimp_update.formula(f, nw~.)
  else nw <- ergm.getnetwork(f)

  m <- ergm.getmodel(f, nw, response=response,...)
  parnames <- param_names(m, canonical=FALSE)
  nparams <- nparam(m, canonical=FALSE)
  coefnames <- param_names(m, canonical=TRUE)
  ncoefs <- nparam(m, canonical=TRUE)

  selection <- unwhich(
    switch(mode(a$which),
           character = match(a$which, parnames),
           logical = which(rep(a$which, length.out=nparams)),
           numeric = a$which),
    nparams)
  if(length(which)) selection[which] <- TRUE
  
  offset.coef <- rep(a$coef, length.out=sum(selection))

  coef0 <- .constrain_init(m, rep(0, nparams))
  coef0[selection] <- offset.coef
    
  Clist <- ergm.Cprepare(nw, m, response=response)

  inputs <- pack.Clist_as_num(Clist)
  
  gs <- ergm.emptynwstats.model(m)
  
  params <- rep(list(NULL), sum(!selection))
  names(params) <- parnames[!selection]
  
  list(name="passthrough_term", coef.names = coefnames, inputs=inputs, dependence=!is.dyad.independent(m), emptynwstats = gs,
       params=params,
       map = function(x, n, ...){
         coef0[!selection] <- x
         ergm.eta(coef0, m$etamap)
       },
       gradient = function(x, n, ...){
         coef0[!selection] <- x
         ergm.etagrad(coef0, m$etamap)[!selection,,drop=FALSE]
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
#' @examples
#' data(sampson)
#' samplike[1,2] <- NA
#' samplike[4.1] <- NA
#' sm <- as.matrix(samplike)
#'
#' tst <- function(x,y){
#'   mapply(identical, x, y)
#' }
#' 
#' stopifnot(all(tst(as.logical(as.matrix(symmetrize(samplike, "weak"))), sm | t(sm))),
#'           all(tst(as.logical(as.matrix(symmetrize(samplike, "strong"))), sm & t(sm))),
#'           all(tst(c(as.matrix(symmetrize(samplike, "upper"))), sm[cbind(c(pmin(row(sm),col(sm))),c(pmax(row(sm),col(sm))))])),
#'           all(tst(c(as.matrix(symmetrize(samplike, "lower"))), sm[cbind(c(pmax(row(sm),col(sm))),c(pmin(row(sm),col(sm))))])))
#' @export
symmetrize.network <- function(x, rule=c("weak","strong","upper","lower"), ...){
  rule <- match.arg(rule)
  el <- rbind(cbind(as.edgelist(x),TRUE),
              cbind(as.edgelist(is.na(x)), NA))
  merge.el <- function(el1, el2){
    # "Encode" NAs as TRUE, TRUE, as FALSE.
    el1[,3] <- is.na(el1[,3])
    el2[,3] <- is.na(el2[,3])
    els <- merge(el1, el2, by.x = c(1:2), by.y = c(1:2), all=TRUE, suffixes = c("th","ht"), sort=FALSE)
    # Now, NA represents FALSE, TRUE NA, and FALSE, TRUE. "Decode".
    els <- within(els,{
      V3th <- ifelse(is.na(V3th), FALSE, ifelse(V3th, NA, TRUE))
      V3ht <- ifelse(is.na(V3ht), FALSE, ifelse(V3ht, NA, TRUE))
    })
    els[els$V1<=els$V2,,drop=FALSE]
  }
  el <- switch(rule,
               weak = {
                 els <- merge.el(el, el[,c(2,1,3)])
                 el <- cbind(els[,1:2,drop=FALSE], els$V3th | els$V3ht)
                 el[is.na(el[,3])|el[,3],,drop=FALSE]
               },                 
               strong = {
                 els <- merge.el(el, el[,c(2,1,3)])
                 el <- cbind(els[,1:2,drop=FALSE], els$V3th & els$V3ht)
                 el[is.na(el[,3])|el[,3],,drop=FALSE]
               },
               upper = el[el[,1]<=el[,2],,drop=FALSE],
               lower = el[el[,1]>=el[,2],,drop=FALSE]
               )
  
  o <- network.initialize(network.size(x), directed=FALSE, bipartite=x%n%"bipartite", loops=has.loops(x), hyper=is.hyper(x), multiple=is.multiplex(x))
  el[,3] <- is.na(el[,3])
  colnames(el) <- c("tails", "heads", "na")
  o <- network.edgelist(el, o, ignore.eval=FALSE, names.eval="na")
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
  
  if(length(f)==2) f <- nonsimp_update.formula(f, nw~.)

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

