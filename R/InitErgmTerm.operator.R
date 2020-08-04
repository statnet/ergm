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
#' @template response
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
wrap.ergm_model <- function(m, nw, response=NULL, namewrap = identity){
  if(!is.null(namewrap)){
    offsettheta <- m$etamap$offsettheta
    offsetmap <- m$etamap$offsetmap
    coef.names <- namewrap(param_names(m, canonical=TRUE))
    ## coef.names[offsetmap] <- paste0("offset(", coef.names[offsetmap], ")")
    minpar <- m$etamap$mintheta
    maxpar <- m$etamap$maxtheta
    # Empty network statistics
    emptynwstats <- summary(m, NULL, response=response)
    if(all(emptynwstats==0)) emptynwstats <- NULL

    # Curved model
    if(is.curved(m)){
      map <- function(x, n, ...){
        ergm.eta(x, m$etamap)
      }
      gradient <- function(x, n, ...){
        ergm.etagrad(x, m$etamap)
      }
      params <- rep(list(NULL), nparam(m))
      names(params) <- namewrap(param_names(m, canonical=FALSE))
      ## names(params)[offsettheta] <- paste0("offset(", names(params)[offsettheta], ")")
    }else map <- gradient <- params <- NULL
  }else{
    minpar <- maxpar <- offsettheta <- offsetmap <- coef.names <- emptynwstats <- map <- gradient <- params <- NULL
  }

  list(map=map, gradient=gradient, params=params, minpar=minpar, maxpar=maxpar, coef.names=coef.names, emptynwstats=emptynwstats, dependence=!is.dyad.independent(m), offsettheta=offsettheta, offsetmap=offsetmap)
}

#' Combine an operator term's and a subterm's name in a standard fashion.
#'
#' @param opname Name of the operator (or an abbreviation thereof).
#' @param opargs A character vector describing arguments passed to the operator (excluding the model); if lengths exceeds one, will be concatenated with commas.
#'
#' @return A function with 1 argument, `subterms`, returning a character vector of length equal to the length of `subterms` wrapped in the operator's name and arguments appropriately.
#' @keywords internal
#' @export
mk_std_op_namewrap <- function(opname, opargs=NULL){
  prefix <-
    if(is.null(opargs)) paste0(opname, "~")
    else paste0(opname, "(", paste(opargs, collapse=","), ")~")
  function(subterms) paste0(prefix, subterms)
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
  
  c(list(name="passthrough_term", submodel=m),
    wrap.ergm_model(m, nw, response, mk_std_op_namewrap('passthrough')))
}

InitErgmTerm.Label <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "label", "pos"),
                      vartypes = c("formula", "character,function,formula,list", "character"),
                      defaultvalues = list(NULL, NULL, "("),
                      required = c(TRUE, TRUE, FALSE))

  f <- a$formula
  if(length(f)==2) f <- nonsimp_update.formula(f, nw~.)
  else nw <- ergm.getnetwork(f)

  m <- ergm_model(f, nw, response=response,...)

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
    wrap.ergm_model(m, nw, response, renamer))
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

  c(list(name="_submodel_term", submodel=m),
    wrap.ergm_model(m, nw, response, NULL))
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

  af <- a$formula
  c(list(name="submodel_test_term", auxiliaries = trim_env(~.submodel(af),"af")),
    wrap.ergm_model(m, nw, response, mk_std_op_namewrap('submodel.test')))
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

  list(name="_summary_term", submodel=m,
       wrap.ergm_model(m, nw, response, NULL))
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

  af <- a$formula
  list(name="summary_test_term", coef.names="summ.test", inputs=c(nparam(m)), auxiliaries=trim_env(~.summary(af),"af"),
       wrap.ergm_model(m, nw, response, NULL))
}

InitErgmTerm.F <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "form"),
                      vartypes = c("formula", "formula"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, TRUE))
  form <- a$form
  f <- a$formula
  if(length(f)==2) f <- nonsimp_update.formula(f, nw~.)
  else nw <- ergm.getnetwork(f)

  m <- ergm_model(f, nw,...)
  
  form.name <- despace(deparse(ult(form)))
  auxiliaries <- trim_env(~.filter.formula.net(form), "form")
  
  c(list(name="on_filter_formula_net",
         submodel = m,
         auxiliaries=auxiliaries),
    wrap.ergm_model(m, nw, response, mk_std_op_namewrap("F", form.name)))
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

  if(!is.dyad.independent(m) || nparam(m)!=1) ergm_Init_abort("The filter test formula must be dyad-independent and have exactly one statistc.")

  nw[,] <- FALSE
  gs <- summary(m, nw, response=response)
  if(gs!=0) ergm_Init_abort("At this time, the filter test term must have the property that its dyadwise components are 0 for 0-valued relations. This limitation may be removed in the future.")
  
  c(list(name="_filter_formula_net", submodel=m),
    wrap.ergm_model(m, nw, response, NULL))
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
  
  offset.coef <- rep(a$coef, length.out=sum(selection))

  coef0 <- .constrain_init(m, rep(0, nparams))
  coef0[selection] <- offset.coef

  params <- rep(list(NULL), sum(!selection))
  names(params) <- parnames[!selection]

  c(list(name="passthrough_term", submodel=m),
    modifyList(wrap.ergm_model(m, nw, response),
               list(coef.names = coefnames,
                    params=params,
                    map = function(x, n, ...){
                      coef0[!selection] <- x
                      ergm.eta(coef0, m$etamap)
                    },
                    gradient = function(x, n, ...){
                      coef0[!selection] <- x
                      ergm.etagrad(coef0, m$etamap)[!selection,,drop=FALSE]
                    }
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

  auxiliaries <- trim_env(~.undir.net(rule), "rule")
  
  c(list(name="on_undir_net",
         submodel = m,
         auxiliaries=auxiliaries),
    modifyList(wrap.ergm_model(m, nw, response, mk_std_op_namewrap('Undir', rule)),
               list(dependence=!is.dyad.independent(m) || rule%in%c("weak","strong"))))
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

  ms <- lapply(fs, ergm_model, nw=nw, response=response, ...)

  curved <- ms[[1]]$etamap$curved
  for(i in seq_len(nf-1L)+1L){
    m <- ms[[i]]
    if(!identical(curved, m$etamap$curved)) ergm_Init_inform("Model ", i, " in the list appears to be curved, and its mapping differs from that of the first model; the first model's mapping will be used.")
  }

  nstats <-  ms %>% map_int(nparam, canonical=TRUE)

  wl <- list()
  for(i in seq_along(fs)){
    f <- fs[[i]]
    w <- if(length(f)==2) 1 else eval_lhs.formula(f)
    if(!is.matrix(w)) w <- diag(w, nstats[i])
    wl[[i]] <- w
  }

  nparams <- wl %>% map_int(nrow)

  if(length(curved) && !all(nparams==nstats[1])) ergm_Init_abort("Specified weights produce different number of output statistics different from those expected by the curved effects in Model 1.")

  if(!all_identical(nparams)) ergm_Init_abort("Specified models and weights appear to differ in lengths of output statistics.")
  nparam <- nparams[1]

  inputs <- unlist(wl%>%map(t))
  inputs <- c(nf, length(inputs), inputs)

  label <- if(length(a$label)==1) paste0(a$label,seq_len(nparam)) else a$label
  coef.names <- mk_std_op_namewrap("Sum")(label)

  wms <- lapply(ms, wrap.ergm_model, nw, response)
  if(is.curved(ms[[1L]])){
    label <- if(length(a$label)==1L) paste0(a$label,seq_along(wms[[1L]]$params)) else attr(a$label,"curved")
    names(wms[[1L]]$params) <- mk_std_op_namewrap("Sum")(label)
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

  if(any(unlist(map(wms, "offsettheta"))) || any(unlist(map(wms, "offsetmap")))) ergm_Init_warn(paste0("Sum operator does not propagate offset() decorators."))
  
  c(list(name="Sum", coef.names = coef.names, inputs=inputs, submodels=ms, emptynwstats=gs, dependence=dependence),
    wms[[1L]][c("map", "gradient", "params", "minpar", "maxpar")])
}

InitErgmTerm.S <- function(nw, arglist, response=NULL, ...){
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
  f <- a$formula
  if(length(f)==3) nw <- ergm.getnetwork(f)
  snw <- nw
  snw[,] <- 0
  snw <- get.inducedSubgraph(snw, tailsel, if(type=="bipartite") headsel)
  if(NVL(snw%n%"bipartite", FALSE)) snw %n% "directed" <- FALSE # Need to do this because snw is a "directed bipartite" network. I hope it doesn't break anything.

  if(length(f)==2) f <- nonsimp_update.formula(f, snw~.)

  m <- ergm_model(f, snw,...)

  auxiliaries <- trim_env(~.subgraph.net(tailsel, headsel), c("tailsel","headsel"))

  selname <- if(tailname==headname) tailname else paste0(tailname,',',headname)

  c(list(name="on_subgraph_net",
         submodel = m,
         auxiliaries=auxiliaries),
    wrap.ergm_model(m, snw, response, mk_std_op_namewrap('S',selname)))
}
