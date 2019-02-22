#' @describeIn to_ergm_Cdouble
#'
#' @return 
#' The `character` method takes a character vector of length 1 and returns a numerical vector encoding it. It concatenates the following:
#' * the number of characters in the input string and
#' * the characters in the input string encoded using ASCII encoding.
#' 
#' This is intended to be decoded by `unpack_str_as_double()` C routines.
#'
#' @export
to_ergm_Cdouble.character <- function(x, ...) c(nchar(x), strtoi(charToRaw(x), 16L))

#' @describeIn to_ergm_Cdouble
#'
#' @return
#' The `ergm_model` method returns a numeric vector concatenating the following:
#' * number of terms in the model;
#' * length of and encoded string of term names;
#' * length of and encoded string of library names; and
#' * vector of intputs to the model.
#' 
#' This is intended to be decoded by `unpack_*Model_as_double()` C routines.
#'
#' @export
to_ergm_Cdouble.ergm_model <- function(x, ...){
  x <- ergm_Clist(x)
  fnames <- to_ergm_Cdouble(x$fnamestring)
  snames <- to_ergm_Cdouble(x$snamestring)
  c(x$nterms, fnames, snames, x$inputs)
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

  m <- ergm_model(f, nw, response=response,...)
  inputs <- to_ergm_Cdouble(m)
  
  gs <- summary(m)
  
  c(list(name="passthrough_term", coef.names = paste0('passthrough(',m$coef.names,')'), inputs=inputs, dependence=!is.dyad.independent(m), emptynwstats = gs),
    passthrough.curved.ergm_model(m, function(x) paste0('passthrough(',x,')')))
}

InitErgmTerm.Label <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "label", "pos"),
                      vartypes = c("formula", "character,function", "character"),
                      defaultvalues = list(NULL, NULL, "("),
                      required = c(TRUE, TRUE, FALSE))

  if(is.character(a$label)){
    pos <- match.arg(a$pos, c("prepend","replace", "(" ,"append"))

    renamer <- switch(pos,
                      prepend=function(x) paste0(a$label,x),
                      replace=function(x) a$label,
                      `(`=function(x) paste0(a$label,"(",x,")"),
                      append=function(x) paste0(x,a$label))
  }else{
    renamer <- a$label
  }

  f <- a$formula
  if(length(f)==2) f <- nonsimp_update.formula(f, nw~.)
  else nw <- ergm.getnetwork(f)

  m <- ergm_model(f, nw, response=response,...)
  inputs <- to_ergm_Cdouble(m)

  gs <- summary(m)

  c(list(name="passthrough_term", coef.names = renamer(m$coef.names), inputs=inputs, dependence=!is.dyad.independent(m), emptynwstats = gs),
    passthrough.curved.ergm_model(m, renamer))
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

  m <- ergm_model(f, nw, response=response,...)
  inputs <- to_ergm_Cdouble(m)

  gs <- summary(m)

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

  m <- ergm_model(f, nw, response=response,...)
  inputs <- to_ergm_Cdouble(m)

  gs <- summary(m)
  
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

  m <- ergm_model(f, nw, response=response,...)
  inputs <- to_ergm_Cdouble(m)

  gs <- summary(m)

  list(name="_summary_term", coef.names = c(), inputs=c(inputs,gs), dependence=!is.dyad.independent(m))
}


## A term to test (very verbosely) the correctness of the .summary() auxiliary.

InitErgmTerm.summary.test <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  f <- a$formula
  if(length(f)==2) f <- nonsimp_update.formula(f, nw~.)
  else nw <- ergm.getnetwork(f)

  m <- ergm_model(f, nw, response=response,...)
  
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

  m <- ergm_model(f, nw,...)
  inputs <- to_ergm_Cdouble(m)
  
  gs <- summary(m)

  form.name <- deparse(ult(form))
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
  
  m <- ergm_model(f, nw, response=response,...)

  if(!is.dyad.independent(m) || nparam(m)!=1) stop("The filter test formula must be dyad-independent and have exactly one statistc.")
  inputs <- to_ergm_Cdouble(m)

  gs <- summary(m)
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

  m <- ergm_model(f, nw, response=response,...)
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
    
  inputs <- to_ergm_Cdouble(m)
  
  gs <- summary(m)
  
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
#'   symmetrized; see [sna::symmetrize()] for details; for the
#'   [`network`] method, it can also be a function or a list; see
#'   Details.
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
#' A method for [`network`] objects, which preserves network and vertex attributes, and handles edge attributes.
#'
#' @details The [`network`] method requires more flexibility, in order
#'   to specify how the edge attributes are handled. Therefore, `rule`
#'   can be one of the following types: \describe{
#' 
#' \item{a character vector}{The string is interpreted as in
#' [sna::symmetrize()]. For edge attributes, `"weak"` takes the
#' maximum value and `"strong"` takes the minimum
#' value" for ordered attributes, and drops the unordered.}
#'
#' \item{a function}{ The function is evaluated on a [`data.frame`]
#' constructed by joining (via [merge()]) the edge [`tibble`] with all
#' attributes and `NA` indicators with itself reversing tail and head
#' columns, and appending original columns with `".th"` and the
#' reversed columns with `".ht"`. It is then evaluated for each
#' attribute in turn, given two arguments: the data frame and the name
#' of the attribute.}
#'
#' \item{a list}{The list must have exactly one unnamed element, and
#' the remaining elements must be named with the names of edge
#' attributes. The elements of the list are interpreted as above,
#' allowing each edge attribute to be handled differently. Unnamed
#' arguments are dropped. }
#' 
#' }
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
  if(!is.directed(x)) return(x)

  TH <- c(".tail",".head")
  HT <- c(".head",".tail")
  
  METACOLS <- c(".tail", ".head", ".eid", "na",".eid.th","na.th",".eid.ht","na.ht")
  
  #' @importFrom tibble as_tibble
  el <- as_tibble(x, attrnames=TRUE, na.rm=FALSE)
  elle <- merge(el, el, by.x = TH, by.y = HT, all=TRUE, suffixes = c(".th",".ht"), sort=FALSE)

  if(!is.list(rule)){
    eattr <- setdiff(list.edge.attributes(x), METACOLS)
    rule <- c(list(rule), rep(list(rule), length(eattr)))
    names(rule) <- c("", eattr)
  }
  
  rule.edges <- which(names(rule)=="")
  if(length(rule.edges)!=1) stop("The list-style argument for rule= must have exactly one unnamed element.")
  rule.edges <- rule[[rule.edges]]


  ## Resolve edges

  NAmap <- function(na) ifelse(is.na(na), FALSE, ifelse(na, NA, TRUE))
  
  elle <-
    if(is.character(rule.edges)){
      keep <-
        switch(match.arg(rule.edges, eval(formals(symmetrize.network)$rule)),
               weak = NAmap(elle$na.th) | NAmap(elle$na.ht),
               strong = NAmap(elle$na.th) & NAmap(elle$na.ht),
               lower = elle$.tail>=elle$.head & NAmap(elle$na.th),
               upper = elle$.tail<elle$.head & NAmap(elle$na.th)
               )
      elle$na <- is.na(keep)
      elle[is.na(keep) | keep,]
    }else{
      rule.edges(elle)
    }

  tmp <- elle$.tail
  elle$.tail <- with(elle, pmin(.tail,.head))
  elle$.head <- with(elle, pmax(tmp,.head))
  elle <- elle[!duplicated(elle[,TH,drop=FALSE]),]
  
  for(attr in setdiff(names(el), METACOLS)){
    r <- NVL(rule[[attr]], rule.edges)
    elle[[attr]] <- ERRVL(try(
      if(is.character(r)){
        th <- paste0(attr,".th")
        ht <- paste0(attr,".ht")
        switch(match.arg(r, eval(formals(symmetrize.network)$rule)),
               pmax =,
               max =, 
               weak = pmax(elle[[th]], elle[[ht]], na.rm=TRUE),
               pmin =,
               min =,
               strong = pmin(elle[[th]], elle[[ht]], na.rm=TRUE),
               lower = elle[[ht]],
               upper = elle[[th]])
      }else{
        r(elle, attr)
      },silent=TRUE), NULL)
  }
  
  o <- network.initialize(network.size(x), directed=FALSE, bipartite=x%n%"bipartite", loops=has.loops(x), hyper=is.hyper(x), multiple=is.multiplex(x))
  add.edges(o, elle$.tail, elle$.head)
  for(attr in c(setdiff(names(el), METACOLS), "na")){
    set.edge.attribute(o, attr, elle[[attr]])
  }
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

  m <- ergm_model(f, nw,...)
  inputs <- to_ergm_Cdouble(m)
  
  gs <- summary(m)

  auxiliaries <- ~.undir.net(rule)
  
  c(list(name="undir",
         coef.names = paste0('Undir(',m$coef.names,',',rule,')'),
         inputs=c(which(RULES==rule),inputs),
         dependence=!is.dyad.independent(m) || rule%in%c("weak","strong"),
         emptynwstats = gs,
         auxiliaries=auxiliaries),
    passthrough.curved.ergm_model(m, function(x) paste0('Undir(',x,',',rule,')')))
}

InitErgmTerm.Sum <- function(nw, arglist, response=NULL,...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formulas", "label"),
                      vartypes = c("list,formula", "character"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, TRUE))

  fs <- a$formulas
  if(is(fs,"formula")) fs <- list(fs)
  nf <- length(fs)
  
  ms <- lapply(fs, function(f){
    m <- ergm_model(f, nw, response=response,...)
    if(is.curved(m)) ergm_Init_inform("Model ", sQuote(deparse(f,500)), " appears to be curved. Its canonical parameters will be used.")
    list(model = m,
         inputs = to_ergm_Cdouble(m),
         gs = summary(m))
  })

  nstats <-  ms %>% map("model") %>% map_int(nparam, canonical=TRUE)

  wl <- list()
  for(i in seq_along(fs)){
    f <- fs[[i]]
    w <- if(length(f)==2) 1 else eval_lhs.formula(f)
    if(!is.matrix(w)) w <- diag(w, nstats[i])
    wl[[i]] <- w
  }

  nparams <- wl %>% map_int(nrow)

  if(!all_identical(nparams)) ergm_Init_abort("Specified models and weights appear to differ in lengths of output statistics.")
  nparam <- nparams[1]

  inputs <- unlist(wl%>%map(t))
  inputs <- c(nf, length(inputs), inputs, ms %>% map("inputs") %>% unlist())

  label <- if(length(a$label)==1) paste0(a$label,seq_len(nparam)) else a$label
  coef.names <- paste0('Sum:',label)
  
  gs <- lst(x = wl,
            y = ms %>% map("gs")) %>%
    pmap(`%*%`) %>%
    reduce(`+`) %>%
    c()
  
  list(name="Sum", coef.names = coef.names, inputs=inputs, dependence=!all(ms %>% map("model") %>% map_lgl(is.dyad.independent)), emptynwstats = gs)
}
