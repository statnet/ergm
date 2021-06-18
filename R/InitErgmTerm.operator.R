#  File R/InitErgmTerm.operator.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################
#' Wrap a submodel's curved, empty network statistics, and extended
#' state (read-only) specification (if present) for output from an
#' `InitErgmTerm` or `InitWtErgmTerm`.
#'
#' Given a `ergm` model and (optionally) a function with which to wrap
#' parameter names, wrap the calls to its `ergm.eta()` and
#' `ergm.etagrad()` into `map()` and `gradient()` functions, similarly
#' with the `params` element; wrap empty network statistics; wrap
#' indicator of dyadic independence; and wrap offset indicators.
#'
#' @param m An `ergm_model` object.
#' @param nw A `network` object.
#' @param namewrap An optional function taking a character vector and
#'   returning a character vector of the same length, called on the
#'   model's canonical and curved parameter names to wrap them. Set to
#'   `NULL` for auxiliary terms to avoid generating elements not
#'   relevant to auxiliaries.
#'
#' @return a list with elements `map`, `gradient`, `params`,
#'   `emptynwstats`, `dependence`, `offsettheta`, and `offsetmap`,
#'   suitable for concatenating with an `InitErgmTerm` or
#'   `InitWtErgmTerm` output list (possibly after modification).
#' @keywords internal
#' @export wrap.ergm_model
wrap.ergm_model <- function(m, nw, namewrap = identity){
  if(!is.null(namewrap)){
    offsettheta <- m$etamap$offsettheta
    coef.names <- namewrap(param_names(m, canonical=TRUE))
    minpar <- m$etamap$mintheta
    maxpar <- m$etamap$maxtheta
    # Empty network statistics
    emptynwstats <- summary(m, NULL)
    if(all(emptynwstats==0)) emptynwstats <- NULL

    # Curved model
    if(is.curved(m)){
      map <- function(x, ...){
        ergm.eta(x, m$etamap)
      }
      gradient <- function(x, ...){
        ergm.etagrad(x, m$etamap)
      }
      params <- rep(list(NULL), nparam(m))
      names(params) <- namewrap(param_names(m, canonical=FALSE))
      ## names(params)[offsettheta] <- paste0("offset(", names(params)[offsettheta], ")")
    }else map <- gradient <- params <- NULL
  }else{
    minpar <- maxpar <- offsettheta <- coef.names <- emptynwstats <- map <- gradient <- params <- NULL
  }

  list(map=map, gradient=gradient, params=params, minpar=minpar, maxpar=maxpar, coef.names=coef.names, emptynwstats=emptynwstats, dependence=!is.dyad.independent(m), offset=offsettheta)
}

#' Combine an operator term's and a subterm's name in a standard fashion.
#'
#' @param opname Name of the operator (or an abbreviation thereof).
#' @param opargs A character vector describing arguments passed to the operator (excluding the model); if lengths exceeds one, will be concatenated with commas.
#'
#' @return A function with 1 required argument, `subterms` and one optional argument, `subargs`, returning a character vector of length equal to the length of `subterms` wrapped in the operator's name and arguments appropriately.
#' @keywords internal
#' @export
ergm_mk_std_op_namewrap <- function(opname, opargs=NULL){
  if(is.null(opargs)) function(subterms, subargs=NULL) NVL3(subargs, paste0(opname, "(", ., ")~", subterms), paste0(opname, "~", subterms))
  else function(subterms, subargs=NULL) NVL3(subargs, paste0(opname, "(", paste0(opargs, collapse=","), ",",., ")~", subterms), paste0(opname, "(", paste(opargs, collapse=","), ")~", subterms))
}

#' Extended states for submodels
#'
#' @description `ergm_propagate_ext.encode()` is a convenience
#'   function to propagate the extended state encoder to submodels if
#'   they have any.
#'
#' @param submodel the [`ergm_model`] to which the encoders should be
#'   propagated.
#'
#' @note `ergm_propagate_ext.encode` should only be used when the
#'   operator term does not modify the network and provides an
#'   `x_function` on the C level that does appropriate propagation and
#'   handles any return values.
#'
#' @return `ergm_propagate_ext.encode` returns a list with one
#'   element, `ext.encode` containing a function that follows the
#'   extended state encoder API and simply returns a list of the
#'   subterms extended state encodings.
#'
#' @examples
#' \dontrun{
#' # Typical usage:
#' InitErgmTerm.ABC <- function(nw, arglist, ...){
#'   [... implementation ...]
#'   m <- ergm_model([... etc. ...])
#'   c(list(name = "abc", inputs=1:3, submodel=m),
#'     ergm_propagate_ext.encode(m),
#'     wrap.ergm_model(nw, m)
#'   )
#' }
#' }
#'
#' @keywords internal
#' @export
ergm_propagate_ext.encode <- function(submodel) {
  has_ext <- !sapply(lapply(submodel$terms, `[[`, "ext.encode"), is.null)

  if (any(has_ext)) list(ext.encode = function(el, nw0)
    lapply(submodel$terms, function(trm) {
      if (!is.null(trm$ext.encode)) trm$ext.encode(el=el, nw0=nw0)
    }))
}

#' @rdname ergm_propagate_ext.encode
#'
#' @description `ergm_no_ext.encode()` checks if a submodel contains
#'   terms that use extended states and stops with an informative
#'   error message if any do.
#'
#' @keywords internal
#' @export
ergm_no_ext.encode <- function(submodel) {
  has_ext <- !sapply(lapply(submodel$terms, `[[`, "ext.encode"), is.null)
  ext_names <- sapply(lapply(submodel$terms[has_ext], `[[`, "call"), deparse, width.cutoff=500)
  if (any(has_ext)) ergm_Init_abort(paste0("This operator term is incompatible with subterms ", paste.and(sQuote(ext_names)), " due to their use of the extended state API. This limitation may be removed in the future."))
}

## Creates a submodel that does exactly what the model terms passed to
## it would have done.
##
InitErgmTerm.Passthrough <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)
  
  c(list(name="passthrough_term", submodel=m),
    ergm_propagate_ext.encode(m),
    wrap.ergm_model(m, nw, ergm_mk_std_op_namewrap('Passthrough')))
}

InitErgmTerm.Label <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "label", "pos"),
                      vartypes = c("formula", "character,function,formula,list", "character"),
                      defaultvalues = list(NULL, NULL, "("),
                      required = c(TRUE, TRUE, FALSE))


  m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)

  if(is.character(a$label)){
    pos <- match.arg(a$pos, c("prepend","replace", "(" ,"append"))
    
    renamer <- switch(pos,
                      prepend=function(x) paste0(a$label,x),
                      replace={
                        if(is.curved(m)){
                          if(!is.list(a$label) || length(a$label)!=2) ergm_Init_abort("For a curved ERGM, replacement label must be a list of length 2, giving the canonical and the curved names, respectively, with NULL to leave alone.")
                          ll <- NVL(a$label[[1]], param_names(m, canonical=TRUE))
                          lc <- NVL(a$label[[1]], param_names(m, canonical=FALSE))
                        }else{
                          ll <- NVL(a$label, param_names(m))
                        }
                        function(x){
                          if(length(x)==nparam(m,canonical=TRUE)) ll
                          else if(length(x)==nparam(m,canonical=FALSE)) lc
                          else ergm_Init_abort("Incorrect length for replacement label vector(s).")
                        }
                      },
                      `(`=function(x) paste0(a$label,"(",x,")"),
                      append=function(x) paste0(x,a$label))
  }else{
    #' @importFrom purrr as_mapper
    renamer <- as_mapper(a$label)
  }

  c(list(name="passthrough_term", submodel=m),
    ergm_propagate_ext.encode(m),
    wrap.ergm_model(m, nw, renamer))
}


## Creates a submodel that tracks the given formula.
InitErgmTerm..submodel <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)
  ergm_no_ext.encode(m)

  c(list(name="_submodel_term", submodel=m),
    wrap.ergm_model(m, nw, NULL))
}

## Tests .submodel().

InitErgmTerm.submodel.test <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)
  ergm_no_ext.encode(m)

  af <- a$formula
  c(list(name="submodel_test_term", auxiliaries = trim_env(~.submodel(af),"af")),
    wrap.ergm_model(m, nw, ergm_mk_std_op_namewrap('submodel.test')))
}

## An auxiliary that exports a double vector that contains the current
## summary statistics of the network for the terms given in the
## formula.

InitErgmTerm..summary <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)

  c(list(name="_summary_term", submodel=m),
    ergm_propagate_ext.encode(m),
    wrap.ergm_model(m, nw, NULL))
}


## A term to test (very verbosely) the correctness of the .summary() auxiliary.

InitErgmTerm.summary.test <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)

  af <- a$formula
  list(name="summary_test_term", coef.names="summ.test", inputs=c(nparam(m)), auxiliaries=trim_env(~.summary(af),"af"),
       wrap.ergm_model(m, nw, NULL))
}

## An auxiliary that exports a model and a double vector that contains
## the current summary statistics of the network for the terms given
## in the formula (or an already initialized model).

InitErgmTerm..submodel_and_summary <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula,ergm_model"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  m <- if(is(a$formula, "formula")) ergm_model(a$formula, nw, ..., offset.decorate=FALSE) else a$formula
  ergm_no_ext.encode(m)

  c(list(name="_submodel_and_summary_term", coef.names = c(), submodel=m),
    wrap.ergm_model(m, nw, NULL))
}

InitErgmTerm.F <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "filter"),
                      vartypes = c("formula", "formula"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, TRUE))

  filter <- a$filter
  m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)
  ergm_no_ext.encode(m)

  filter.name <- despace(deparse(ult(filter)))
  auxiliaries <- trim_env(~.filter.formula.net(filter), "filter")

  c(list(name="on_filter_formula_net",
         submodel = m,
         auxiliaries=auxiliaries),
    wrap.ergm_model(m, nw, ergm_mk_std_op_namewrap("F", filter.name)))
}

InitErgmTerm..filter.formula.net <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)
  ergm_no_ext.encode(m)

  if(!is.dyad.independent(m) || nparam(m)!=1) ergm_Init_abort("The filter test formula must be dyad-independent and have exactly one statistc.")

  nw[,] <- FALSE
  gs <- summary(m, nw)
  if(gs!=0) ergm_Init_abort("At this time, the filter test term must have the property that its dyadwise components are 0 for 0-valued relations. This limitation may be removed in the future.")
  
  c(list(name="_filter_formula_net", submodel=m),
    wrap.ergm_model(m, nw, NULL))
}

InitErgmTerm.Offset <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "coef", "which"),
                      vartypes = c("formula", "numeric", "logical,numeric,character"),
                      defaultvalues = list(NULL, 0, TRUE),
                      required = c(TRUE, FALSE, FALSE))

  m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)

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
  
  offset.coef <- rep(a$coef, length.out=sum(selection))

  coef0 <- rep(NA, nparams)
  coef0[selection] <- offset.coef

  params <- rep(list(NULL), sum(!selection))
  names(params) <- parnames[!selection]

  c(list(name="passthrough_term", submodel=m),
    ergm_propagate_ext.encode(m),
    replace(wrap.ergm_model(m, nw),
            c("coef.names", "params", "map", "gradient", "offset"),
            list(coef.names = coefnames,
                    params=params,
                    map = function(x, n, ...){
                      coef0[!selection] <- x
                      ergm.eta(coef0, m$etamap)
                    },
                    gradient = function(x, n, ...){
                      coef0[!selection] <- x
                      ergm.etagrad(coef0, m$etamap)[!selection,,drop=FALSE]
                    },
                    offset = m$etamap$offsettheta[!selection]
                    )))
}

#' Return a symmetrized version of a binary network
#'
#' @param x an object representing a network.
#' @param rule a string specifying how the network is to be
#'   symmetrized; see [sna::symmetrize()] for details; for the
#'   [`network`] method, it can also be a function or a list; see
#'   Details.
#' @param ... additional arguments to [sna::symmetrize()].
#'
#' @note This was originally exported as a generic to overwrite
#'   [sna::symmetrize()]. By developer's request, it has been renamed;
#'   eventually, `sna` or `network` packages will export the generic
#'   instead.
#' @export
ergm_symmetrize <- function(x, rule=c("weak","strong","upper","lower"), ...){
  UseMethod("ergm_symmetrize")
}

#' @describeIn ergm_symmetrize
#'
#' The default method, passing the input on to [sna::symmetrize()].
#' 
#' @export
ergm_symmetrize.default <- function(x, rule=c("weak","strong","upper","lower"), ...){
  sna::symmetrize(x, rule=rule, ...)
}

#' @describeIn ergm_symmetrize
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
#' stopifnot(all(tst(as.logical(as.matrix(ergm_symmetrize(samplike, "weak"))), sm | t(sm))),
#'           all(tst(as.logical(as.matrix(ergm_symmetrize(samplike, "strong"))), sm & t(sm))),
#'           all(tst(c(as.matrix(ergm_symmetrize(samplike, "upper"))),
#'                   sm[cbind(c(pmin(row(sm),col(sm))),c(pmax(row(sm),col(sm))))])),
#'           all(tst(c(as.matrix(ergm_symmetrize(samplike, "lower"))),
#'                   sm[cbind(c(pmax(row(sm),col(sm))),c(pmin(row(sm),col(sm))))])))
#' @export
ergm_symmetrize.network <- function(x, rule=c("weak","strong","upper","lower"), ...){
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
        switch(match.arg(rule.edges, eval(formals(ergm_symmetrize.network)$rule)),
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
        switch(match.arg(r, eval(formals(ergm_symmetrize.network)$rule)),
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

InitErgmTerm.Symmetrize <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist, directed = TRUE,
                      varnames = c("formula", "rule"),
                      vartypes = c("formula", "character"),
                      defaultvalues = list(NULL, "weak"),
                      required = c(TRUE, FALSE))
  RULES <- c("weak","strong","upper","lower")
  rule <- match.arg(a$rule, RULES)

  if(is.directed(nw)) nw <- ergm_symmetrize(nw, rule)
  m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)
  ergm_no_ext.encode(m)

  auxiliaries <- trim_env(~.undir.net(rule), "rule")
  
  c(list(name="on_undir_net",
         submodel = m,
         auxiliaries=auxiliaries),
    modifyList(wrap.ergm_model(m, nw, ergm_mk_std_op_namewrap('Symmetrize', rule)),
               list(dependence=!is.dyad.independent(m) || rule%in%c("weak","strong"))))
}

InitErgmTerm.Sum <- function(nw, arglist,...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formulas", "label"),
                      vartypes = c("list,formula", "character,function"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, TRUE))

  fs <- a$formulas
  if(is(fs,"formula")) fs <- list(fs)
  nf <- length(fs)

  ms <- lapply(fs, ergm_model, nw=nw, ..., offset.decorate=FALSE)

  curved <- ms[[1]]$etamap$curved
  for(i in seq_len(nf-1L)+1L){
    m <- ms[[i]]
    if(!identical(curved, m$etamap$curved)) ergm_Init_inform("Model ", i, " in the list appears to be curved, and its mapping differs from that of the first model; the first model's mapping will be used.")
  }

  nstats <-  ms %>% map_int(nparam, canonical=TRUE)

  wl <- list()
  for(i in seq_along(fs)){
    f <- fs[[i]]
    w <- if(length(f)==2) 1
         else if(!is.character(lhs <- eval_lhs.formula(f))) lhs
         else switch(match.arg(lhs, c("sum","mean")),
                     sum = matrix(1, 1, nstats[i]),
                     mean = matrix(1/nstats[i], 1, nstats[i]))
    if(!is.matrix(w)) w <- diag(w, nstats[i])
    wl[[i]] <- w
  }

  nparams <- wl %>% map_int(nrow)

  if(length(curved) && !all(nparams==nstats[1])) ergm_Init_abort("Specified weights produce different number of output statistics different from those expected by the curved effects in Model 1.")

  if(!all_identical(nparams)) ergm_Init_abort("Specified models and weights appear to differ in lengths of output statistics.")
  nparam <- nparams[1]

  inputs <- unlist(wl%>%map(t))
  inputs <- c(nf, length(inputs), inputs)


  if(is.function(a$label)){
    cns <- lapply(ms, param_names, canonical=TRUE)
    if(length(cns) == 1) cns <- cns[[1]]
    cn <- a$label(cns)
  }else cn <- a$label
  cn.asis <- inherits(cn, "AsIs")

  cn <- if(length(cn)==1L && nparam>1L) paste0(cn, seq_len(nparam)) else cn
  if(length(cn) != nparam) ergm_Init_abort(paste0(sQuote("label="), " argument for statistics has or results in length ", length(cn), ", should be ", nparam, "."))
  coef.names <- if(cn.asis) cn else ergm_mk_std_op_namewrap("Sum")(cn)

  wms <- lapply(ms, wrap.ergm_model, nw)
  if(is.curved(ms[[1L]])){
    ncparam <- length(wms[[1L]]$params)

    if(is.function(a$label)){
      pns <- lapply(ms, param_names, canonical=FALSE)
      if(length(pns) == 1) pns <- pns[[1]]
      pn <- a$label(pns)
    }else pn <- NVL(attr(a$label,"curved"), a$label)
    pn.asis <- inherits(pn, "AsIs")

    pn <- if(length(pn)==1L && ncparam>1L) paste0(pn, seq_len(ncparam)) else pn
    if(length(pn) != ncparam) ergm_Init_abort(paste0(sQuote("label="), " argument for curved parameters has or results in length ", length(pn), ", should be ", ncparam, "."))
    names(wms[[1L]]$params) <- if(pn.asis) pn else ergm_mk_std_op_namewrap("Sum")(pn)
  }

  gss <- map(wms, "emptynwstats")
  gs <-
    if(all(map_lgl(gss, is.null))){ # Linear combination of 0s is 0.
      NULL
    }else{ # All numeric or NULL
      gs0 <- map_lgl(gss, is.null)
      lst(x = wl[!gs0],
          y = gss[!gs0]) %>%
        pmap(`%*%`) %>%
        reduce(`+`) %>%
        c()
    }

  dependence <- any(map_lgl(wms, "dependence"))

  if(any(unlist(map(wms, "offset")))){
    if(length(wl) == 1 && diff(w <- dim(wl[[1]])) == 0 && # There is only one formula, its weights are a square matrix,
       !any(w[!diag(ncol(w))])) # and all its off-diagonal elements are 0,
      offset <- wms[[1L]]$offset # then offsets are safe to propagate.
    else{
      offset <- FALSE
      ergm_Init_warn(paste0("Sum operator does not propagate offset() decorators unless there is only one formula and its statistics are simply scaled."))
    }
  }else offset <- FALSE
  
  c(list(name="Sum", coef.names = coef.names, inputs=inputs, submodels=ms, emptynwstats=gs,
         dependence=dependence, offset=offset,
         ext.encode = if(ms %>% map("terms") %>% unlist(FALSE) %>% map("ext.encode") %>% compact %>% length)
                        function(el, nw0)
                          lapply(ms, function(submodel)
                            lapply(submodel$terms, function(trm){
                              if(!is.null(trm$ext.encode)) trm$ext.encode(el=el, nw0=nw0)
                            }))),
    wms[[1L]][c("map", "gradient", "params", "minpar", "maxpar")])
}

InitErgmTerm.S <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "attrs"),
                      vartypes = c("formula", ERGM_VATTR_SPEC),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, TRUE))

  ### Obtain tail and head selectors.
  spec <- a$attrs
  if(!is(spec, "formula")) spec <- call("~",spec) # Embed into formula.
  if(length(spec)==2) spec <- call("~", spec[[2]], spec[[2]]) # Duplicate if one-sided.

  class(spec) <- "formula"
  environment(spec) <- environment(a$attrs)

  # Obtain subformulas for tail and head.
  tailspec <- spec[-3]
  headspec <- spec[-2]

  # Obtain the boolean indicators or numeric indices. If the network
  # is bipartite in the first place, expect bipartite indices.
  bip <- as.integer(NVL(nw %n% "bipartite", 0L))

  tailsel <- ergm_get_vattr(tailspec, nw, accept="index", bip="b1")
  tailname <- attr(tailsel, "name")

  headsel <- ergm_get_vattr(headspec, nw, accept="index", bip="b2")
  headname <- attr(headsel, "name")

  # Convert to numeric selectors.
  if(is.logical(tailsel)) tailsel <- which(tailsel)
  if(is.logical(headsel)) headsel <- which(headsel) + bip

  tailsel <- as.integer(tailsel)
  headsel <- as.integer(headsel)
  
  # TODO: Check if 1-node graphs cause crashes.
  if(length(tailsel)==0 || length(headsel)==0) ergm_Init_abort("Empty subgraph selected.")

  type <- if(is.directed(nw)) "directed" else "undirected"
  if(bip){
    if(max(tailsel)>bip || min(headsel)<=bip)
      ergm_Init_abort("Invalid vertex subsets selected for a bipartite graph.")
    type <- "bipartite"
  }else{
    if(!identical(tailsel,headsel)){ # Rectangular selection: output bipartite.
      if(length(intersect(tailsel,headsel))) ergm_Init_abort("Vertex subsets constructing a bipartite subgraph must have disjoint ranges.")
      type <- "bipartite"
    }
  }

  ### Construct an empty network with the correct structure.
  snw <- nw
  snw[,] <- 0
  snw <- get.inducedSubgraph(snw, tailsel, if(type=="bipartite") headsel)
  if(NVL(snw%n%"bipartite", FALSE)) snw %n% "directed" <- FALSE # Need to do this because snw is a "directed bipartite" network. I hope it doesn't break anything.

  m <- ergm_model(a$formula, snw, ..., offset.decorate=FALSE)
  ergm_no_ext.encode(m)

  auxiliaries <- trim_env(~.subgraph.net(tailsel, headsel), c("tailsel","headsel"))

  selname <- if(tailname==headname) tailname else paste0(tailname,',',headname)

  c(list(name="on_subgraph_net",
         submodel = m,
         auxiliaries=auxiliaries),
    wrap.ergm_model(m, snw, ergm_mk_std_op_namewrap('S',selname)))
}


InitErgmTerm.Parametrise <- InitErgmTerm.Parametrize <- InitErgmTerm.Curve <- function(nw, arglist,...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "params", "map", "gradient", "minpar", "maxpar", "cov"),
                      vartypes = c("formula", "character,list", "function,numeric,character", "function,matrix,character", "numeric", "numeric", ""),
                      defaultvalues = list(NULL, NULL, NULL, NULL, -Inf, +Inf, NULL),
                      required = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE))

  m <- ergm_model(a$formula, nw=nw, ..., offset.decorate=FALSE)
  p <- nparam(m, canonical=FALSE)

  params <- if(is.character(a$params)) setNames(rep(list(NULL), length(a$params)), a$params) else a$params
  q <- length(params)

  minpar <- rep(a$minpar, length.out=q)
  maxpar <- rep(a$maxpar, length.out=q)

  emap <- # map() is used by purrr.
    if(is.numeric(a$map)){
      NVL(a$gradient) <- "linear" # TODO: Replace with a constructed matrix.
      function(...) a$map
    }else if(is.character(a$map)){
      switch(match.arg(a$map, c("rep")),
             rep = {
               NVL(a$gradient) <- "linear" # TODO: Replace with a constructed matrix.
               function(x, n, ...) rep_len(x, n)
             })
    }else a$map

  if(is.null(a$gradient)) ergm_Init_abort(paste0("The ", sQuote("gradient"), " argument must be supplied unless ", sQuote("map"), " is of a special type."))

  gradient <-
    if(is.matrix(a$gradient)) function(...) a$gradient
    else if(is.character(a$gradient))
      switch(match.arg(a$gradient, c("linear")),
             linear = { # Gradient is constant.

               x <- (deInf(minpar) + deInf(maxpar))/2
               u <- pmin(x+1/2, maxpar)
               l <- pmax(x-1/2, minpar)
               d <- u-l

               G <- vapply(seq_len(q), function(i){
                 (emap(ifelse(i==seq_len(q), u, x), p, a$cov)
                   - emap(ifelse(i==seq_len(q), l, x), p, a$cov))/d[i]
               }, numeric(p)) %>% t
               function(...) G
             })
    else a$gradient

  # Make sure the output dimensions are correct.
  test.param <- (deInf(minpar) + deInf(maxpar))/2
  test.map <- emap(test.param, p, a$cov)
  if(length(test.map)!=p) ergm_Init_abort(paste0("Model expects ", p, " parameters, but the map function returned a vector of length ", length(test.map), "."))
  test.gradient <- gradient(test.param, p, a$cov)
  if(!identical(dim(test.gradient),c(q,p))) ergm_Init_abort(paste0("Mapping of ", q, " to ", p, " parameters expected, but the gradient function returned an object with dimension (", paste0(dim(test.gradient), collapse=","), ")."))

  wm <- wrap.ergm_model(m, nw)

  if(any(unlist(map(wm, "offsettheta"))) || any(unlist(map(wm, "offsetmap")))) ergm_Init_warn(paste0("Curve operator does not propagate offset() decorators."))

  if(is.curved(m)){
    wm$map <- function(x, n, ...) wm$map(emap(x, n, ...)) # Composition.
    wm$gradient <- function(x, n, ...) gradient(x, n, ...) %*% wm$gradient(emap(x, n, ...)) # Chain rule.
  }else{
    wm$map <- emap
    wm$gradient <- gradient
  }

  wm$cov <- a$cov
  wm$params <- params
  wm$minpar <- minpar
  wm$maxpar <- maxpar
  wm$offset <- logical(q)

  c(list(name="passthrough_term", submodel=m), wm,
    ergm_propagate_ext.encode(m))
}

InitErgmTerm.Log <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "log0"),
                      vartypes = c("formula", "numeric"),
                      defaultvalues = list(NULL, -1/sqrt(.Machine$double.eps)),
                      required = c(TRUE, FALSE))

  m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)
  ergm_no_ext.encode(m)
  log0 <- rep_len(a$log0, nparam(m, canonical=TRUE))

  wm <- wrap.ergm_model(m, nw, ergm_mk_std_op_namewrap('Log'))
  NVL(wm$emptynwstats) <- numeric(nparam(m, canonical=TRUE))
  wm$emptynwstats <- ifelse(wm$emptynwstats==0, log0, log(wm$emptynwstats))

  c(list(name="Log", inputs=log0, auxiliaries=~.submodel_and_summary(m)), wm)
}

InitErgmTerm.Exp <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)
  ergm_no_ext.encode(m)

  wm <- wrap.ergm_model(m, nw, ergm_mk_std_op_namewrap('Exp'))
  wm$emptynwstats <- wm$emptynwstats %>% NVL(numeric(nparam(m, canonical=TRUE))) %>% exp

  c(list(name="Exp", auxiliaries=~.submodel_and_summary(m)), wm)
}

InitErgmTerm.Prod <- function(nw, arglist, ..., env=baseenv()){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formulas", "label"),
                      vartypes = c("list,formula", "character,function"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, TRUE))
  formulas <- if(is(a$formulas, "formula")) list(a$formulas) else a$formulas
  formulas <- lapply(formulas, function(f) nonsimp_update.formula(f, .~Log(~.)))
  formulas <- lapply(formulas, function(f) {
    if(length(f)==3 && is.character(f[[2]]))
      f[[2]] <- switch(match.arg(f[[2]], c("prod","geomean")),
                       prod = "sum",
                       geomean = "mean")
    f
  })
  a$formulas <- if(length(formulas)==1) formulas[[1]] else formulas
  cl <- call("Exp", as.formula(call("~",call("Sum", a$formulas, a$label)),env=env))
  call.ErgmTerm(cl, env, nw, ...)
}
