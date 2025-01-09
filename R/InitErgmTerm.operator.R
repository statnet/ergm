#  File R/InitErgmTerm.operator.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
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
#' `namewrap` also controls how dyadic dependence flag is propagated
#' for auxiliaries. If `NULL`, it is propagated; if not, the
#' auxiliaries are ignored and only terms's dyadic dependence is
#' propagated.
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

    dependence <- !is.dyad.independent(m, ignore_aux=TRUE)
  }else{
    minpar <- maxpar <- offsettheta <- coef.names <- emptynwstats <- map <- gradient <- params <- NULL
    dependence <- !is.dyad.independent(m, ignore_aux=FALSE)
  }

  list(map=map, gradient=gradient, params=params, minpar=minpar, maxpar=maxpar, coef.names=coef.names, emptynwstats=emptynwstats, dependence=dependence, offset=offsettheta)
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
  if (any(has_ext)) ergm_Init_stop("This operator term is incompatible with subterms ", paste.and(sQuote(ext_names)), " due to their use of the extended state API. This limitation may be removed in the future.")
}

## Creates a submodel that does exactly what the model terms passed to
## it would have done.
##
InitErgmTerm.Passthrough <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "label", "submodel"),
                      vartypes = c("formula", "logical", "logical"),
                      defaultvalues = list(NULL, FALSE, TRUE),
                      required = c(TRUE, FALSE, FALSE))

  if(a$submodel){
    m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)

    c(list(name="passthrough_term", submodel=m),
      ergm_propagate_ext.encode(m),
      wrap.ergm_model(m, nw, if(a$label) ergm_mk_std_op_namewrap('Passthrough') else identity))
  }else{
    m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE, terms.only=TRUE)

    namewrap <- if(a$label) ergm_mk_std_op_namewrap('Passthrough') else identity
    param_names(m) <- list(namewrap(param_names(m, canonical = FALSE)), namewrap(param_names(m, canonical = TRUE)))
    m
  }
}

#' @templateVar name Label
#' @title Modify terms' coefficient names
#' @description This operator evaluates `formula` without modification, but modifies its coefficient and/or parameter names based on `label` and `pos` .
#'
#' @usage
#' # binary: Label(formula, label, pos)
#' @template ergmTerm-formula
#' @param label a character vector specifying the label for the terms, a [`list`] of two character vectors (see Details), or a function through which term names are mapped (or a [`as_mapper`] -style formula).
#' @param pos controls how `label` modifies the term names: one of `"prepend"` , `"replace"` , `"append"` , or `"("` , with the latter wrapping the term names in parentheses like a function call with name specified by `label` .
#'
#' @details If `pos == "replace"`:
#'
#' * Elements for which `is.na(label) == TRUE` are preserved.
#'
#' * If the model is curved, `label=` can be a either function/mapper
#'   or a [`list`] with two elements, the first element giving the
#'   curved (model) parameter names and second giving the canonical
#'   parameter names. `NULL` leaves the respective name unchanged.
#'
#' @template ergmTerm-general
#'
#' @concept operator
InitErgmTerm.Label <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "label", "pos"),
                      vartypes = c("formula", "character,function,formula,list", "character"),
                      defaultvalues = list(NULL, NULL, "("),
                      required = c(TRUE, TRUE, FALSE))

  m <- ergm_model(a$formula, nw, ..., offset.decorate = FALSE, terms.only = TRUE)
  cu <- param_names(m, canonical = FALSE)
  ca <- param_names(m, canonical = TRUE)

  if(is.character(a$label) || is.list(a$label)){
    pos <- match.arg(a$pos, c("prepend","replace", "(" ,"append"))
    
    new <- switch(pos,
                  prepend = list(paste0(a$label, cu), paste0(a$label, ca)),
                  replace =
                    if(is.curved(m)){
                      if(!is.list(a$label) || length(a$label)!=2) ergm_Init_stop("For a curved ERGM, replacement label must be a list of length 2, giving the curved and the canonical names, respectively, with NULL to leave alone.")
                      list(NVL(a$label[[1]], NA), NVL(a$label[[2]], NA))
                    }else rep(list(NVL(a$label, cu)), 2),
                  `(` = list(paste0(a$label,"(",cu,")"), paste0(a$label,"(",ca,")")),
                  append = list(paste0(cu, a$label), paste0(ca, a$label))
                  )
  }else{
    #' @importFrom purrr as_mapper
    renamer <- as_mapper(a$label)
    new <- list(cu,ca) %>% map(renamer)
  }

  param_names(m) <- new
  m
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

# Roxygen converts the name to FALSE-ergmTerm if the name is not specified. This generates a warning but is probably the best we can do

#' @templateVar name 'F
#' @title Filtering on arbitrary one-term model
#' @description Evaluates the given `formula` on a network constructed by
#'   taking \eqn{y} and removing any edges for which
#'   \eqn{f_{i,j}(y_{i,j}) = 0}{f[i,j] (y[i,j])=0} .
#'
#' @usage
#' # binary: F(formula, filter)
#' @template ergmTerm-formula
#' @param filter must contain one binary `ergm` term, with
#'   the following properties:
#'   - dyadic independence;
#'   - dyadwise contribution of 0 for a 0-valued dyad.
#'
#'   Formally, this means that it is expressable as
#'   \deqn{g(y) = \sum_{i,j} f_{i,j}(y_{i,j}),}{sum[i,j]
#'   f[i,j] (y[i,j]),} where for all \eqn{i}, \eqn{j}, and \eqn{y},
#'   \eqn{f_{i,j}(y_{i,j})} for which \eqn{f_{i,j}(0)=0}{f[i,j] (0)=0}.
#'   For convenience, the term in specified can be a part of a simple logical or comparison operation: (e.g., `~!nodematch("A")` or `~abs("X")>3`),
#'   which filters on \eqn{f_{i,j}(y_{i,j}) \bigcirc 0}{f[i,j] (y[i,j]) \%OP\% 0} instead.
#'
#' @template ergmTerm-general
#'
#' @concept operator
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

  OPS <- c("!", "==", "!=", ">" , "<", ">=" , "<=")
  UNARY <- c("!")
  if(is.call(ult(a$formula)) && (op <- as.character(ult(a$formula)[[1]])) %in% OPS){
    iinputs <- match(op, OPS)
    inputs <- if(! op%in%UNARY) eval(ult(a$formula)[[3]], environment(a$formula))
    ult(a$formula) <- ult(a$formula)[[2]]
  }else{
    iinputs <- 0L
    inputs <- NULL
  }

  m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)
  ergm_no_ext.encode(m)

  if(!is.dyad.independent(m) || nparam(m)!=1) ergm_Init_stop("The filter test formula must be dyad-independent and have exactly one statistic.")

  nw[,] <- FALSE
  gs <- summary(m, nw)
  if(gs!=0) ergm_Init_stop("At this time, the filter test term must have the property that its dyadwise components are 0 for 0-valued relations. This limitation may be removed in the future.")
  
  c(list(name="_filter_formula_net", submodel=m, iinputs=iinputs, inputs=inputs),
    wrap.ergm_model(m, nw, NULL))
}

#' @templateVar name Offset
#' @title Terms with fixed coefficients
#' @description This operator is analogous to the `offset()` wrapper, but the
#'   coefficients are specified within the term and the curved ERGM
#'   mechanism is used internally.
#'
#' @usage
#' # binary: Offset(formula, coef, which)
#' @template ergmTerm-formula
#' @param coef coefficients to the formula
#' @param which used to specify which of the parameters in the formula are fixed. It can be a logical vector (recycled as needed), a numeric vector of indices of parameters to be fixed, or a character vector of parameter names.
#'
#' @template ergmTerm-general
#'
#' @concept operator
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
#'   eventually, \CRANpkg{sna} or `network` packages will export the generic
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

#' @templateVar name Symmetrize
#' @title Evaluation on symmetrized (undirected) network
#' @description Evaluates the terms in `formula` on an undirected network
#'   constructed by symmetrizing the LHS network using one of four rules:
#'   
#'   1. "weak" A tie \eqn{(i,j)} is present in the constructed
#'   network if the LHS network has either tie \eqn{(i,j)} or
#'   \eqn{(j,i)} (or both).
#'   2. "strong" A tie \eqn{(i,j)} is present in the constructed
#'   network if the LHS network has both tie \eqn{(i,j)} and tie
#'   \eqn{(j,i)} .
#'   3. "upper" A tie \eqn{(i,j)} is present in the constructed
#'   network if the LHS network has tie \eqn{(\min(i,j),\max(i,j))} :
#'   the upper triangle of the LHS network.
#'   4. "lower" A tie \eqn{(i,j)} is present in the constructed
#'   network if the LHS network has tie \eqn{(\max(i,j),\min(i,j))} :
#'   the lower triangle of the LHS network.
#'
#' @usage
#' # binary: Symmetrize(formula, rule="weak")
#' @template ergmTerm-formula
#' @param rule one of `"weak"`, `"strong"`, `"upper"`, `"lower"`
#'
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept operator
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

#' @templateVar name Sum
#' @title A sum (or an arbitrary linear combination) of one or more formulas
#' @description This operator sums up the RHS statistics of the input formulas elementwise.
#'   
#'
#' @details Note that each formula must either produce the same number of
#'   statistics or be mapped through a matrix to produce the same
#'   number of statistics.
#'   
#'   A single formula is also permitted. This can be useful if one
#'   wishes to, say, scale or sum up the statistics returned by a formula.
#'
#'   Offsets are ignored unless there is only one formula and the transformation only scales the statistics (i.e., the effective transformation matrix is diagonal).
#'   
#'   Curved models are supported, subject to some limitations. In particular, the first model's etamap will be used, overwriting the others. If `label` is not of length 1, it should have an `attr` -style attribute `"curved"` specifying the names for the curved parameters.
#'
#' @usage
#' # binary: Sum(formulas, label)
#' @param formulas a list (constructed using [list()] or [c()]) of [ergm()]-style formulas whose RHS gives the statistics to be evaluated, or a single formula.
#'
#'   If a formula in the list has an LHS, it is interpreted as follows:
#'   - a numeric scalar: Network statistics of this formula will be multiplied by this.
#'   - a numeric vector: Corresponding network statistics of this formula will be multiplied by this.
#'   - a numeric matrix: Vector of network statistics will be pre-multiplied by this.
#'   - a character string: One of several predefined linear combinations. Currently supported presets are as follows:
#'     - `"sum"` Network statistics of this formula will be summed up; equivalent to `matrix(1,1,p)` , where `p` is the length of the network statistic vector.
#'     - `"mean"` Network statistics of this formula will be averaged; equivalent to `matrix(1/p,1,p)` , where `p` is the length of the network statistic vector.
#' @param label used to specify the names of the elements of the resulting term sum vector. If `label` is a character vector of length 1,
#'   it will be recycled with indices appended. If a function is specified, `formulas` parameter names are extracted and their list of character vectors is passed `label`.
#"   (For convenience, if only one formula is given, just a character vector is passed. Lastly, if `label` or result of its function call is an [`AsIs`] object, it is not wrapped in `Sum~...`.)
#'
#' @template ergmTerm-general
#'
#' @concept operator
InitErgmTerm.Sum <- function(nw, arglist,...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formulas", "label"),
                      vartypes = c("list,formula", "character,function,AsIs"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, TRUE))

  fs <- a$formulas
  if(is(fs,"formula")) fs <- list(fs)
  nf <- length(fs)

  ms <- lapply(fs, ergm_model, nw=nw, ..., offset.decorate=FALSE)

  curved <- ms[[1]]$etamap$curved
  for(i in seq_len(nf-1L)+1L){
    m <- ms[[i]]
    if(!identical(curved, m$etamap$curved)) ergm_Init_message("Model ", i, " in the list appears to be curved, and its mapping differs from that of the first model; the first model's mapping will be used.")
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

  if(length(curved) && !all(nparams==nstats[1])) ergm_Init_stop("Specified weights produce different number of output statistics different from those expected by the curved effects in Model 1.")

  if(!all_identical(nparams)) ergm_Init_stop("Specified models and weights appear to differ in lengths of output statistics.")
  nparam <- nparams[1]

  inputs <- unlist(wl%>%map(t))

  if(is.function(a$label)){
    cns <- lapply(ms, param_names, canonical=TRUE)
    if(length(cns) == 1) cns <- cns[[1]]
    cn <- a$label(cns)
  }else cn <- a$label
  cn.asis <- inherits(cn, "AsIs")

  cn <- if(length(cn)==1L && nparam>1L) paste0(cn, seq_len(nparam)) else cn
  if(length(cn) != nparam) ergm_Init_stop(sQuote("label="), " argument for statistics has or results in length ", length(cn), ", should be ", nparam, ".")
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
    if(length(pn) != ncparam) ergm_Init_stop(sQuote("label="), " argument for curved parameters has or results in length ", length(pn), ", should be ", ncparam, ".")
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
      ergm_Init_warning("Sum operator does not propagate offset() decorators unless there is only one formula and its statistics are simply scaled.")
    }
  }else offset <- FALSE
  
  c(list(name="Sum", coef.names = coef.names, inputs=inputs, iinputs=nf, submodels=ms, emptynwstats=gs,
         dependence=dependence, offset=offset,
         ext.encode = if(ms %>% map("terms") %>% unlist(FALSE) %>% map("ext.encode") %>% compact %>% length)
                        function(el, nw0)
                          lapply(ms, function(submodel)
                            lapply(submodel$terms, function(trm){
                              if(!is.null(trm$ext.encode)) trm$ext.encode(el=el, nw0=nw0)
                            }))),
    wms[[1L]][c("map", "gradient", "params", "minpar", "maxpar")])
}

#' @templateVar name S
#' @title Evaluation on an induced subgraph
#' @description This operator takes a two-sided forumla `attrs` whose LHS gives the attribute or attribute function for which tails and heads will be used to construct the induced subgraph. They must evaluate either to a logical vector equal in length to the number of tails (for LHS) and heads (for RHS) indicating which nodes are to be used to induce the subgraph or a numeric vector giving their indices. 
#'
#' @details As with indexing vectors, the logical vector will be recycled to the size of the network or the size of the appropriate bipartition, and negative indices will deselect vertices.
#'   
#'   When the two sets are identical, the induced subgraph retains the directedness of the original graph. Otherwise, an undirected bipartite graph is induced.
#'
#' @usage
#' # binary: S(formula, attrs)
#' @template ergmTerm-formula
#' @param attrs a two-sided formula to be used. A one-sided formula (e.g., `~A` ) is symmetrized (e.g., `A~A` ).
#'
#' @template ergmTerm-general
#'
#' @concept operator
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
  if(length(tailsel)==0 || length(headsel)==0) ergm_Init_stop("Empty subgraph selected.")

  type <- if(is.directed(nw)) "directed" else "undirected"
  if(bip){
    if(max(tailsel)>bip || min(headsel)<=bip)
      ergm_Init_stop("Invalid vertex subsets selected for a bipartite graph.")
    type <- "bipartite"
  }else{
    if(!identical(tailsel,headsel)){ # Rectangular selection: output bipartite.
      if(length(intersect(tailsel,headsel))) ergm_Init_stop("Vertex subsets constructing a bipartite subgraph must have disjoint ranges.")
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


#' @templateVar name Curve
#' @title Impose a curved structure on term parameters
#' @description Arguments may have the same forms as in the API, but for convenience, alternative forms are accepted.
#'   
#'   If the model in `formula` is curved, then the outputs of this operator term's `map` argument will be used as inputs to the curved terms of the `formula` model.
#'
#'   `Curve` is an obsolete alias and may be deprecated and removed in a future release.
#'
#' @usage
#' # binary: Curve(formula, params, map, gradient=NULL, minpar=-Inf, maxpar=+Inf, cov=NULL)
#' @template ergmTerm-formula
#' @param params a named list whose names are the curved parameter names, may also be a character vector with names.
#' @param map the mapping from curved to canonical. May have the following forms:
#' - a `function(x, n, ...)` treated as in the API: called with `x` set to the curved parameter vector, `n` to the length of output expected, and `cov` , if present, passed in `...` . The function must return a numeric vector of length `n` .
#' - a numeric vector to fix the output coefficients, like in an offset.
#' - a character string to select (partially-matched) one of predefined forms. Currently, the defined forms include:
#'   - `"rep"` recycle the input vector to the length of the output vector as a `rep` function would.
#' @param gradient its gradient function. It is optional if `map` is constant or one of the predefined forms; otherwise it must have one of the following forms:
#' - a `function(x, n, ...)` treated as in the API: called with `x` set to the curved parameter vector, `n` to the length of output expected, and `cov` , if present, passed in `...` . The function must return a numeric matrix with `length(params)` rows and `n` columns.
#' - a numeric matrix to fix the gradient; this is useful when map is linear.
#' - a character string to select (partially-matched) one of predefined forms. Currently, the defined forms include:
#'   - `"linear"` calculate the (constant) gradient matrix using finite differences. Note that this will be done only once at the initialization stage, so use only if you are certain `map` is, in fact, linear.
#' @param minpar,maxpar the minimum and maximum allowed curved parameter values. The parameters will be recycled to the appropriate length.
#' @param cov optional
#'
#' @template ergmTerm-general
#'
#' @concept operator
InitErgmTerm.Curve <- function(nw, arglist,...){
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

  if(is.null(a$gradient)) ergm_Init_stop("The ", sQuote("gradient"), " argument must be supplied unless ", sQuote("map"), " is of a special type.")

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
  if(length(test.map)!=p) ergm_Init_stop("Model expects ", p, " parameters, but the map function returned a vector of length ", length(test.map), ".")
  test.gradient <- gradient(test.param, p, a$cov)
  if(!identical(dim(test.gradient),c(q,p))) ergm_Init_stop("Mapping of ", q, " to ", p, " parameters expected, but the gradient function returned an object with dimension (", paste0(dim(test.gradient), collapse=","), ").")

  wm <- wrap.ergm_model(m, nw)

  if(any(unlist(map(wm, "offsettheta"))) || any(unlist(map(wm, "offsetmap")))) ergm_Init_warning("Curve operator does not propagate ", sQuote("offset()"), " decorators.")

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

#' @templateVar name Curve
#' @template ergmTerm-rdname
#' @aliases Parametrise-ergmTerm
#' @usage
#' # binary: Parametrise(formula, params, map, gradient=NULL, minpar=-Inf, maxpar=+Inf,
#' #           cov=NULL)
InitErgmTerm.Parametrise <- InitErgmTerm.Curve

#' @templateVar name Curve
#' @template ergmTerm-rdname
#' @aliases Parametrize-ergmTerm
#' @usage
#' # binary: Parametrize(formula, params, map, gradient=NULL, minpar=-Inf, maxpar=+Inf,
#' #           cov=NULL)
InitErgmTerm.Parametrize <- InitErgmTerm.Curve

#' @templateVar name Log
#' @title Take a natural logarithm of a network's statistic
#' @description Evaluate the terms specified in `formula` and takes a natural (base \eqn{e} ) logarithm of them. Since an ERGM statistic must be finite, `log0` specifies the value to be substituted for `log(0)` . The default value seems reasonable for most purposes.
#'
#' @usage
#' # binary: Log(formula, log0=-1/sqrt(.Machine$double.eps))
#' @template ergmTerm-formula
#' @param log0 the value to be substituted for `log(0)`
#'
#' @template ergmTerm-general
#'
#' @concept operator
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

#' @templateVar name Exp
#' @title Exponentiate a network's statistic
#' @description Evaluate the terms specified in `formula` and exponentiates them with base \eqn{e} .
#'
#' @usage
#' # binary: Exp(formula)
#' @template ergmTerm-formula
#'
#' @template ergmTerm-general
#'
#' @concept operator
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

#' @templateVar name Prod
#' @title A product (or an arbitrary power combination) of one or more formulas
#' @description This operator evaluates a list of formulas whose corresponnding RHS
#'   statistics will be multiplied elementwise. They are required to be nonnegative.
#'   
#' @details    
#' Note that each formula must either produce the same number of
#'   statistics or be mapped through a matrix to produce the same
#'   number of statistics.
#'   
#'   A single formula is also permitted. This can be useful if one
#'   wishes to, say, scale or multiply together the statistics returned by a formula.
#'   
#'   Offsets are ignored unless there is only one formula and the transformation only scales the statistics (i.e., the effective transformation matrix is diagonal).
#'
#'   Curved models are supported, subject to some limitations. In particular, the first model's etamap will be used, overwriting the others. If `label` is not of length 1, it should have an `attr` -style attribute `"curved"` specifying the names for the curved parameters.
#'   
#' @usage
#' # binary: Prod(formulas, label)
#' @param formulas a list (constructed using [list()] or [c()]) of [ergm()]-style formulas whose RHS gives the statistics to be evaluated, or a single formula.
#'
#' If a formula in the list has an LHS, it is interpreted as follows:
#'   - a numeric scalar: Network statistics of this formula will be exponentiated by this.
#'   - a numeric vector: Corresponding network statistics of this formula will be exponentiated by this.
#'   - a numeric matrix: Vector of network statistics will be exponentiated by this using the same pattern as matrix multiplication.
#'   - a character string: One of several predefined multiplicative combinations. Currently supported presets are as follows:
#'     - `"prod"`: Network statistics of this formula will be multiplied together; equivalent to `matrix(1,1,p)` , where `p` is the length of the network statistic vector.
#'     - `"geomean"`: Network statistics of this formula will be geometrically averaged; equivalent to `matrix(1/p,1,p)` , where `p` is the length of the network statistic vector.

#' @param label used to specify the names of the elements of the resulting term product vector. If `label` is a character vector of length 1,
#'   it will be recycled with indices appended. If a function is specified, `formulas` parameter names are extracted and their list of character vectors is passed `label`.
#"   (For convenience, if only one formula is given, just a character vector is passed. Lastly, if `label` or result of its function call is an [`AsIs`] object, it is not wrapped in `Prod~...`.
#'
#' @template ergmTerm-general
#'
#' @note The current implementation piggybacks on the `Log` , `Exp` , and `Sum` operators, essentially `Exp(~Sum(~Log(formula), label))` . This may result in loss of precision, particularly for extremely large or small statistics. The implementation may change in the future.
#'
#' @concept operator
InitErgmTerm.Prod <- function(nw, arglist, ..., env=baseenv()){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formulas", "label"),
                      vartypes = c("list,formula", "character,function,AsIs"),
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
  trm <- call.ErgmTerm(cl, env, nw, ...)
  renamer <-
    if(inherits(a$label, "AsIs")) function(x) sub("^Exp~", "", x)
    else function(x) sub("^Exp~Sum~", "Prod~", x)

  trm$coef.names %<>% renamer()
  if(!is.null(trm$params)) names(trm$params) %<>% renamer()

  trm
}


#' @templateVar name For
#' @title A [`for`] operator for terms
#' @description This operator evaluates the formula given to it,
#'   substituting the specified loop counter variable with each
#'   element in a sequence.
#'
#' @usage
#' # binary: For(...)
#' @param ... in any order, \itemize{
#'
#'   \item one *unnamed* one-sided [ergm()]-style formula with the
#'   terms to be evaluated, containing one or more placeholders
#'   \var{VAR} *and*
#'
#'   \item one or more *named* expressions of the form \code{\var{VAR}
#'   = \var{SEQ}} specifying the placeholder and its range. See
#'   Details below.
#'
#' }
#'
#' @details Placeholders are specified in the style of
#'   `foreach::foreach()`, as \code{\var{VAR} = \var{SEQ}}. \var{VAR}
#'   can be any valid \R variable name, and \var{SEQ} can be a vector,
#'   a list, a function of one argument, or a one-sided formula.  The
#'   vector or list will be used directly, whereas a function will be
#'   called with the network as its argument to produce the list, and
#'   the formula will be used analogously to [purrr::as_mapper()], its
#'   RHS evaluated in an environment in which the network itself will
#'   be accessible as `.` or `.nw`.
#'
#'   If more than one named expression is given, they will be expanded
#'   as one would expect in a nested [`for`] loop: earlier expressions
#'   will form the outer loops and later expressions the inner loops.
#'
#' @template ergmTerm-general
#'
#' @examples
#' #
#' # The following are equivalent ways to compute differential
#' # homophily.
#' #
#'
#' data(sampson)
#' (groups <- sort(unique(samplike%v%"group"))) # Sorted list of groups.
#'
#' # The "normal" way:
#' summary(samplike ~ nodematch("group", diff=TRUE))
#'
#' # One element at a time, specifying a list:
#' summary(samplike ~ For(~nodematch("group", levels=., diff=TRUE),
#'                        . = groups))
#'
#' # One element at a time, specifying a function that returns a list:
#' summary(samplike ~ For(~nodematch("group", levels=., diff=TRUE),
#'                        . = function(nw) sort(unique(nw%v%"group"))))
#'
#' # One element at a time, specifying a formula whose RHS expression
#' # returns a list:
#' summary(samplike ~ For(~nodematch("group", levels=., diff=TRUE),
#'                        . = ~sort(unique(.%v%"group"))))
#'
#' #
#' # Multiple iterators are possible, in any order. Here, absdiff() is
#' # being computed for each combination of attribute and power.
#' #
#'
#' data(florentine)
#'
#' # The "normal" way:
#' summary(flomarriage ~ absdiff("wealth", pow=1) + absdiff("priorates", pow=1) +
#'                       absdiff("wealth", pow=2) + absdiff("priorates", pow=2) +
#'                       absdiff("wealth", pow=3) + absdiff("priorates", pow=3))
#'
#' # With a loop; note that the attribute (a) is being iterated within
#' # power (.):
#' summary(flomarriage ~ For(. = 1:3, a = c("wealth", "priorates"), ~absdiff(a, pow=.)))
#'
#' @concept operator
InitErgmTerm.For <- function(nw, arglist, ...){
  counters <- names(arglist)
  if(length(i <- which(counters=="")) != 1) ergm_Init_stop("Exactly one argument (the model formula) must be unnamed.")
  if(length(loops <- arglist[-i]) < 1) ergm_Init_stop("At least one counter must be provided.")

  loops <- map(loops,
               function(l){
                 if(is(l, "formula")){
                   vlist <- lst(`.`=nw, `.nw`=nw)
                   e <- ult(l)
                   ergm_Init_try({eval(e, envir=vlist, enclos=environment(l))})
                 }else if(is.function(l)){
                   ergm_Init_try({l(nw)})
                 }else l
               })

  valid <- loops %>% map_lgl(is.vector)
  if(!all(valid)) ergm_Init_stop("Loop variable(s) ", paste.and(sQuote(names(loops)[!valid])), " does not contain a valid sequence.")

  terms <- list_rhs.formula(arglist[[i]])

  for(i in rev(seq_along(loops))){
    counter <- names(loops)[i]
    loop <- loops[[i]]
    terms <- do.call(c,
                     map(loop, function(.x){
                       for(j in seq_along(terms)){
                         terms[[j]] <- do.call(substitute, list(terms[[j]], structure(list(.x), names=counter)))
                       }
                       terms
                     }))
  }

  ergm_model(terms, nw, ..., terms.only=TRUE)
}
