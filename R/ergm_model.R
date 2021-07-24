#  File R/ergm_model.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################
#===================================================================================
# This file contains the following 2 functions for creating the 'ergm_model' object
#             <ergm_model>
#             <updatemodel.ErgmTerm>
#===================================================================================

#' Internal representation of an `ergm` network model
#' 
#' These methods are generally not called directly by users, but may
#' be employed by other depending packages.
#' `ergm_model` constructs it from a formula. Each term is
#' initialized via the \code{InitErgmTerm} functions to create a
#' \code{ergm_model} object.
#' @note This API is not to be considered fixed and may change between versions. However, an effort will be made to ensure that the methods of this class remain stable.
#' @param formula An [ergm()]
#' formula of the form \code{network ~ model.term(s)} or \code{~
#' model.term(s)}.
#' @param nw The network of interest; if passed, the LHS of `formula` is ignored. This is the recommended usage.
#' @param silent logical, whether to print the warning messages from the
#' initialization of each model term.
#' @param \dots additional parameters for model formulation
#' @template term_options
#' @param extra.aux a list of auxiliary request formulas required elsewhere; if named, the resulting `slots.extra.aux` will also be named.
#' @param env a throwaway argument needed to prevent conflicts with some usages of `ergm_model`. The initialization environment is *always* taken from the `formula`.
#' @param offset.decorate logical; whether offset coefficient and parameter names should be enclosed in `"offset()"`. 
#' @param object An `ergm_model` object.
#' @return `ergm_model` returns an  `ergm_model` object as a list
#' containing:
#' \item{coef.names}{a vector of coefficient names (deprecated: use [param_names( , canonical=TRUE)] instead)}
#' \item{terms}{a list of terms and 'term components' initialized by the
#' appropriate \code{InitErgmTerm.X} function.}
#' \item{etamap}{the theta -> eta mapping as a list returned from
#' <ergm.etamap>}
#' @seealso [summary.ergm_model()]
#' @keywords internal
#' @export
ergm_model <- function(formula, nw=NULL, silent=FALSE, ..., term.options=list(), extra.aux=list(), env=globalenv(), offset.decorate=TRUE){
  if (!is(formula, "formula"))
    stop("Invalid model formula of class ",sQuote(class(formula)),".", call.=FALSE)
  
  #' @importFrom statnet.common eval_lhs.formula
  if(is.null(nw)) nw <- eval_lhs.formula(formula)

  nw <- ensure_network(nw)
  nw <- as.network(nw, populate=FALSE) # In case it's an ergm_state.
  
  #' @importFrom statnet.common list_rhs.formula
  v<-list_rhs.formula(formula)
  
  formula.env<-environment(formula)
  
  model <- structure(list(coef.names = character(),
                          terms = list(), networkstats.0 = numeric()),
                 class = "ergm_model")

  for (i in 1:length(v)) {
    term <- v[[i]]

    if (is.call(term) && term[[1L]] == "offset"){ # Offset term
      offset <-
        if(length(term)==3) eval(term[[3]], formula.env)
        else TRUE
      term <- term[[2L]]
    }else offset <- FALSE

    ## term is now a call or a name that is not "offset".
    
    if(!is.call(term) && term==".") next
    
    outlist <- call.ErgmTerm(term, formula.env, nw, term.options=term.options, ...)
    
    # If initialization fails without error (e.g., all statistics have been dropped), continue.
    if(is.null(outlist)){
      if(!silent) message("Note: Term ", deparse(v[[i]])," skipped because it contributes no statistics.")
      model$term.skipped <- c(model$term.skipped, TRUE)
      next
    }else model$term.skipped <- c(model$term.skipped, FALSE)

    # Now it is necessary to add the output to the model formula
    model <- updatemodel.ErgmTerm(model, outlist, offset=offset, offset.decorate=offset.decorate, silent=silent)
  }

  model <- ergm.auxstorage(model, nw, term.options=term.options, ..., extra.aux=extra.aux)
  
  model$etamap <- ergm.etamap(model)

  if(offset.decorate){
    if(length(model$etamap$offsetmap)){
      model$coef.names <- ifelse(model$etamap$offsetmap, paste0("offset(",model$coef.names,")"), model$coef.names)

      ol <- split(model$etamap$offsetmap, factor(rep.int(seq_along(model$terms), nparam(model, byterm=TRUE, canonical=TRUE)), levels=seq_along(model$terms)))
      for(i in seq_along(model$terms)){
        pn <- model$terms[[i]]$coef.names
        if(!is.null(pn)) model$terms[[i]]$coef.names <- ifelse(ol[[i]], paste0("offset(",pn,")"), pn)
      }
    }
    if(length(model$etamap$offsettheta)){
      ol <- split(model$etamap$offsettheta, factor(rep.int(seq_along(model$terms), nparam(model, byterm=TRUE, canonical=FALSE)), levels=seq_along(model$terms)))
      for(i in seq_along(model$terms)){
        pn <- names(model$terms[[i]]$params)
        if(!is.null(pn)) names(model$terms[[i]]$params) <- ifelse(ol[[i]], paste0("offset(",pn,")"), pn)
      }
    }
  }

  # I.e., construct a vector of package names associated with the model terms.
  # Note that soname is not the same, since it's not guaranteed to be a loadable package.
  ergm.MCMC.packagenames(unlist(sapply(model$terms, "[[", "pkgname")))

  class(model) <- "ergm_model"
  model
}

#' Locate and call an ERGM term initialization function.
#'
#' A helper function that searches attached and loaded packages for a
#' term with a specifies name, calls it with the specified arguments,
#' and returns the result.
#'
#' @param A term from an [ergm()] formula: typically a [`name`] or a
#'   [`call`].
#' @param env Environment in which it is to be evaluated.
#' @param nw A [`network`] object.
#' @param term.options A list of optional settings such as calculation
#'   tuning options to be passed to the `InitErgmTerm` functions.
#' @param ... Additional term options.
#'
#' @return The list returned by the the `InitErgmTerm` or
#'   `InitWtErgmTerm` function, with package name autodetected if
#'   neede.
#'
#' @keywords internal
#' @export call.ErgmTerm
call.ErgmTerm <- function(term, env, nw, ..., term.options=list()){
  termroot<-if(!is.valued(nw)) "InitErgm" else "InitWtErgm"
  
  if(is.call(term)) { # This term has some arguments; save them.
    args <- term
    args[[1]] <- as.name("list")
  }else args <- list()
  
  termFun <- locate_prefixed_function(term, paste0(termroot,"Term"), "ERGM term", env=env)
  termCall <- termCall(termFun, term, nw, term.options, ..., env=env)

  #Call the InitErgm function in the environment where the formula was created
  # so that it will have access to any parameters of the ergm terms
  out <- eval(termCall, env)
  if(is.null(out)) return(NULL)
  # If SO package name not specified explicitly, autodetect.
  if(is.null(out$pkgname)) out$pkgname <- environmentName(environment(eval(termFun)))
  # Store the term call in the term list (with the term able to override)
  NVL(out$call) <- term
  # If the term requests auxiliaries or is an auxiliary itself,
  # reserve space in the input vector.
  attr(out, "aux.slots") <-
    integer(NVL3(out$auxiliaries, length(list_rhs.formula(.)), 0) + # requests auxiliaries
            (length(out$coef.names)==0)) # is an auxiliary

  # Ensure input vectors are of the correct storage mode. (There is no
  # checking on C level at this time.) Note that as.double() and
  # as.integer() will strip attributes such as ParamBeforeCov and so
  # should not be used.
  storage.mode(out$inputs) <- "double"
  storage.mode(out$iinputs) <- "integer"
  if(!is.null(out$emptynwstats)) storage.mode(out$emptynwstats) <- "double"

  out
}

#' Updates an existing model object to include an initialized `ergm` term
#'
#' @param model the pre-existing model, as created by [`ergm_model`]
#' @param outlist the list describing new term, as returned by `InitErgmTerm.*()`
#' @param offset a vector indicating which parameters in the term should be offsets
#' @param silent whether to suppress messages
#'
#' @return The updated model (with the obvious changes seen below) if
#'   `outlist!=NULL`, else the original model. (Note that this return
#'   is necessary, since terms may be eliminated by giving only 0
#'   statistics, and consequently returning a NULL `outlist`.)
#' @noRd
updatemodel.ErgmTerm <- function(model, outlist, offset=FALSE, offset.decorate=TRUE, silent=FALSE) {
  if (!is.null(outlist)) { # Allow for no change if outlist==NULL
    # Update global model properties.
    nstats <- length(outlist$coef.names)
    npars <- NVL3(outlist$params, length(.), nstats)

    if(!silent && npars == 0 && ((is.numeric(offset) && length(offset)) || (is.logical(offset) && any(offset))))
      message(sQuote("offset()"), " decorator used on term ",sQuote(deparse(outlist$call)), " with no free parameters is meaningless and will be ignored.")
    if(is.numeric(offset)) offset <- unwhich(offset, npars)
    outlist$offset <- offset <- rep(offset, length.out=npars) | NVL(outlist$offset,FALSE)

    model$coef.names <- c(model$coef.names, outlist$coef.names)
    model$minval <- c(model$minval,
                      rep(NVL(outlist$minval, -Inf),
                          length.out=nstats))
    model$maxval <- c(model$maxval,
                      rep(NVL(outlist$maxval, +Inf),
                          length.out=nstats))
    model$terms[[length(model$terms)+1L]] <- outlist
  }
  model
}

#' @describeIn ergm_model A method for concatenating terms of two or more initialized models.
#' @export
c.ergm_model <- function(...){
  l <- list(...)
  o <- l[[1]]

  # Generally useful:
  naux <- function(trm) NVL3(trm$auxiliaries, length(list_rhs.formula(.)), 0L)

  for(m in l[-1]){
    if(is.null(m)) next

    oaux <- o$terms %>% map("coef.names") %>% map_int(length)==0
    maux <- m$terms %>% map("coef.names") %>% map_int(length)==0

    # Sort auxiliaries into sinks (that don't depend on any others) and non-sinks:
    ## TODO: Implement an iterative or graph algorithm that checks *all* auxiliaries for redundancy.
    pos <- 0L
    onaux <- o$terms[oaux] %>% map_int(naux)
    omap <- integer(length(o$terms[oaux]))
    for(i in seq_along(o$terms[oaux]))
      if(!onaux[i]){ # sink
        omap[i] <- attr(o$terms[oaux][[i]],"aux.slots")[1L] <- NA
      }else{ # non-sink
        omap[i] <- attr(o$terms[oaux][[i]],"aux.slots")[1L] <- pos
        pos <- pos+1L
      }

    mnaux <- m$terms[maux] %>% map_int(naux)
    mmap <- integer(length(m$terms[maux]))
    for(i in seq_along(m$terms[maux]))
      if(!mnaux[i]){ # sink
        mmap[i] <- attr(m$terms[maux][[i]],"aux.slots")[1L] <- NA
      }else{ # non-sink
        mmap[i] <- attr(m$terms[maux][[i]],"aux.slots")[1L] <- pos
        pos <- pos+1L
      }

    # Remove redundant auxiliaries, but only among sinks:
    sinks <- unique_aux_terms(c(o$terms[oaux][!onaux], m$terms[maux][!mnaux]))
    
    # To which term in the new auxilary list does each of the current auxiliary pointers to sinks correspond? Note that they are numbered from 0.
    omap[!onaux] <- match_aux_terms(o$terms[oaux][!onaux], sinks) + pos - 1L
    mmap[!mnaux] <- match_aux_terms(m$terms[maux][!mnaux], sinks) + pos - 1L
    
    # Now, give them indices again:
    for(i in seq_along(sinks))
      attr(sinks[[i]],"aux.slots")[[1L]] <- pos + i - 1L

    # Remap client term pointers:
    for(i in seq_along(o$terms[!oaux]))
      if((tnaux <- naux(o$terms[!oaux][[i]]))!=0L)
        attr(o$terms[!oaux][[i]],"aux.slots")[seq_len(tnaux)] <- omap[attr(o$terms[!oaux][[i]],"aux.slots")[seq_len(tnaux)]+1L]
    o$slots.extra.aux <-  map(o$slots.extra.aux, ~omap[.+1L])
    
    for(i in seq_along(m$terms[!maux]))
      if((tnaux <- naux(m$terms[!maux][[i]]))!=0L)
        attr(m$terms[!maux][[i]],"aux.slots")[seq_len(tnaux)] <- mmap[attr(m$terms[!maux][[i]],"aux.slots")[seq_len(tnaux)]+1L]
    m$slots.extra.aux <- map(m$slots.extra.aux, ~mmap[.+1L])

    # Similarly for auxiliaries as clients.
    for(i in seq_along(o$terms[oaux])[onaux!=0L])
      attr(o$terms[oaux][[i]],"aux.slots")[1L+seq_len(onaux[i])] <- omap[attr(o$terms[oaux][[i]],"aux.slots")[1L+seq_len(onaux[i])]+1L]
    for(i in seq_along(m$terms[maux])[mnaux!=0L])
      attr(m$terms[maux][[i]],"aux.slots")[1L+seq_len(mnaux[i])] <- mmap[attr(m$terms[maux][[i]],"aux.slots")[1L+seq_len(mnaux[i])]+1L]

    # New auxiliary list is o-nonsinks, m-nonsinks, and sinks.
    o$terms <- c(o$terms[!oaux], m$terms[!maux], o$terms[oaux][onaux!=0L], m$terms[maux][mnaux!=0L], sinks)

    for(name in c("coef.names",
                  "minval",
                  "maxval",
                  "offset",
                  "term.skipped",
                  "slots.extra.aux"))
      o[[name]] <- c(o[[name]], m[[name]])
  }

  # Check that terms and auxiliaries are properly positioned.
  assert_aux_dependencies(o$terms)

  o$etamap <- ergm.etamap(o)
  o
}

#' @rdname ergm_model
#' @param x object to be converted to an `ergm_model`.
#' @export
as.ergm_model <- function(x, ...)
  UseMethod("as.ergm_model")

#' @rdname ergm_model
#' @export
as.ergm_model.ergm_model <- function(x, ...) x

#' @rdname ergm_model
#' @export
as.ergm_model.formula <- function(x, ...)
  ergm_model(formula=x, ...)

#' @noRd
#' @export
as.ergm_model.NULL <- function(x, ...) NULL


#' Construct a standard term call
#'
#' The call will comprise a function, the network, an argument list
#' (whose composition depends on the API), and a list of additional
#' options, passed either through `term.options` or in `...`.
#'
#' @param f the function to be called.
#' @param args a list of arguments to be passed as a part of the
#'   argument list or a [call()]: if the latter, the first element is
#'   replaced by `list`.
#' @param nw a network.
#' @template term_options
#' @param ... additional options.
#' @noRd
termCall <- function(f, args, nw, term.options, ...){
  if(is.call(args)) args[[1]] <- as.name("list") # This term has some arguments; save them.
  else if(!is.list(args)) args <- list()

  termCall <- as.call(list(f, nw, args))

  #' @importFrom utils modifyList
  dotdotdot <- as.list(getOption("ergm.term")) %>% modifyList(as.list(term.options)) %>% modifyList(list(...))
  # FIXME: There's probably a more elegant way to do this, but this
  # works for both calls and lists.
  for(i in seq_along(dotdotdot)) {
    termCall[[length(termCall)+1]] <- dotdotdot[[i]]
    names(termCall)[length(termCall)] <- names(dotdotdot)[i]
  }

  termCall
}
