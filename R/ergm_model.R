#  File R/ergm_model.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
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
#' @template response
#' @param silent logical, whether to print the warning messages from the
#' initialization of each model term.
#' @param role A hint about how the model will be used. Used primarily for
#' dynamic network models.
#' @param \dots additional parameters for model formulation
#' @param term.options a list of optional settings such as calculation tuning options to be passed to the `InitErgmTerm` functions.
#' @param extra.aux a list of auxiliary request formulas required elsewhere.
#' @param object An `ergm_model` object.
#' @return `ergm_model` returns an  `ergm_model` object as a list
#' containing:
#' \item{ formula}{the formula inputted to
#' \code{\link{ergm_model}}}
#' \item{coef.names}{a vector of coefficient names}
#' \item{offset}{a logical vector of whether each term was "offset", i.e.
#' fixed}
#' \item{terms}{a list of terms and 'term components' initialized by the
#' appropriate \code{InitErgmTerm.X} function.}
#' \item{network.stats0}{NULL always??}
#' \item{etamap}{the theta -> eta mapping as a list returned from
#' <ergm.etamap>}
#' @seealso [summary.ergm_model()]
#' @export
ergm_model <- function(formula, nw=NULL, response=NULL, silent=FALSE, role="static",...,term.options=list(),extra.aux=list()){
  if (!is(formula, "formula"))
    stop("Invalid model formula of class ",sQuote(class(formula)),".", call.=FALSE)
  
  #' @importFrom statnet.common eval_lhs.formula
  if(is.null(nw)) nw <- eval_lhs.formula(formula)

  #' @importFrom utils modifyList
  term.options <- modifyList(as.list(getOption("ergm.term")), as.list(term.options))
  
  #' @importFrom statnet.common list_rhs.formula
  v<-list_rhs.formula(formula)
  
  formula.env<-environment(formula)
  
  model <- structure(list(formula=formula, coef.names = NULL,
                      offset = NULL,
                      terms = NULL, networkstats.0 = NULL, etamap = NULL),
                 class = "ergm_model")
  
  for (i in 1:length(v)) {
    term <- v[[i]]
    
    if (is.call(term) && term[[1]] == "offset"){ # Offset term
      term <- term[[2]]
      model$offset <- c(model$offset,TRUE)
    }else{
      model$offset <- c(model$offset,FALSE)
    }
    ## term is now a call or a name that is not "offset".
    
    if(!is.call(term) && term==".") next
    
    outlist <- call.ErgmTerm(term, formula.env, nw, response=response, role=role, term.options=term.options, ...)
    
    # If initialization fails without error (e.g., all statistics have been dropped), continue.
    if(is.null(outlist)){
      if(!silent) message("Note: Term ", deparse(v[[i]])," skipped because it contributes no statistics.")
      model$term.skipped <- c(model$term.skipped, TRUE)
      next
    }else model$term.skipped <- c(model$term.skipped, FALSE)
    # If the term is an offset, rename the coefficient names and parameter names
    if(model$offset[length(model$offset)]){
      outlist$coef.names <- paste0("offset(",outlist$coef.names,")")
      if(!is.null(outlist$params))
        names(outlist$params) <- paste0("offset(",outlist$params,")")
    }
    # Now it is necessary to add the output to the model formula
    model <- updatemodel.ErgmTerm(model, outlist)
  }

  model <- ergm.auxstorage(model, nw, response=response, term.options=term.options, ..., extra.aux=extra.aux)
  
  model$etamap <- ergm.etamap(model)

  # I.e., construct a vector of package names associated with the model terms.
  # Note that soname is not the same, since it's not guaranteed to be a loadable package.
  ergm.MCMC.packagenames(unlist(sapply(model$terms, "[[", "pkgname")))
  ergm.MCMC.packagenames(unlist(sapply(model$model.aux$terms, "[[", "pkgname")))
  
  class(model) <- "ergm_model"
  model
}


call.ErgmTerm <- function(term, env, nw, response=NULL, role="static", ..., term.options=list()){
  termroot<-if(is.null(response)) "InitErgm" else "InitWtErgm"
  
  if(is.call(term)) { # This term has some arguments; save them.
    args <- term
    args[[1]] <- as.name("list")
  }else args <- list()
  
  termFun<-locate.InitFunction(term, paste0(termroot,"Term"), "ERGM term", env=env)
  termCall<-as.call(list(termFun, nw, args))
  
  dotdotdot <- c(if(!is.null(response)) list(response=response), list(...), list(role=role), term.options)
  for(j in seq_along(dotdotdot)) {
    if(is.null(dotdotdot[[j]])) next
    termCall[[3+j]] <- dotdotdot[[j]]
    names(termCall)[3+j] <- names(dotdotdot)[j]
  }
  #Call the InitErgm function in the environment where the formula was created
  # so that it will have access to any parameters of the ergm terms
  out <- eval(termCall,env)
  # If SO package name not specified explicitly, autodetect.
  if(!is.null(out) && is.null(out$pkgname)) out$pkgname <- environmentName(environment(eval(termFun)))
  out
}

#' @describeIn ergm-deprecated Use `ergm_model` instead.
#' @export ergm.getmodel
ergm.getmodel <- function(object, ...){
  .Deprecated("ergm_model")
  ergm_model(object, ...)
}

#######################################################################
# The <updatemodel.ErgmTerm> function updates an existing model object
# to include an initialized ergm term, X;
#
# --PARAMETERS--
#   model  : the pre-existing model, as created by <ergm_model>
#   outlist: the list describing term X, as returned by <InitErgmTerm.X>
#
# --RETURNED--
#   model: the updated model (with the obvious changes seen below) if
#            'outlist'!=NULL, else
#          the original model; (note that this return is necessary,
#            since terms may be eliminated by giving only 0 statistics,
#            and consequently returning a NULL 'outlist')
#
#######################################################################

updatemodel.ErgmTerm <- function(model, outlist) { 
  if (!is.null(outlist)) { # Allow for no change if outlist==NULL
    model$coef.names <- c(model$coef.names, outlist$coef.names)
    termnumber <- 1+length(model$terms)
    tmp <- attr(outlist$inputs, "ParamsBeforeCov")
    # If the term requests auxiliaries or is an auxiliary itself,
    # reserve space in the input vector. Note that these go before
    # the parameters.
    aux.space <-
      if(!is.null(outlist$auxiliaries)) # requests auxiliaries
        length(list_rhs.formula(outlist$auxiliaries))
      else if(length(outlist$coef.names)==0) 1 # is an auxiliary
      else 0
    outlist$inputs <- c(ifelse(is.null(tmp), 0, tmp)+aux.space,
                        length(outlist$coef.names), 
                        length(outlist$inputs)+aux.space, rep(NA, aux.space), outlist$inputs)
    model$minval <- c(model$minval,
                      rep(NVL(outlist$minval, -Inf),
                          length.out=length(outlist$coef.names)))
    model$maxval <- c(model$maxval,
                      rep(NVL(outlist$maxval, +Inf),
                          length.out=length(outlist$coef.names)))
    model$duration <- max(model$duration,
                          NVL(outlist$duration, FALSE))
    model$terms[[termnumber]] <- outlist
  }
  model
}

#' @describeIn ergm_model A method for concatenating terms of two or more initialized models.
#' @export
c.ergm_model <- function(...){
  l <- list(...)
  o <- l[[1]]
  for(m in l[-1]){
    if(is.null(m)) next

    # Reset auxiliaries' pointers to own slot to NA:
    for(i in seq_along(o$model.aux$terms))
      o$model.aux$terms[[i]]$inputs[4] <- NA
    for(i in seq_along(m$model.aux$terms))
      m$model.aux$terms[[i]]$inputs[4] <- NA

    # Remove redundant auxiliaries:
    uniq.aux.outlists <- unique(c(o$model.aux$terms, m$model.aux$terms))
    
    # To which term in the new auxilary list does each of the current auxiliary pointers correspond? Note that they are numbered from 0.
    omap <- match(o$model.aux$terms, uniq.aux.outlists)-1
    mmap <- match(m$model.aux$terms, uniq.aux.outlists)-1
    
    # Now, give them indices again:
    for(i in seq_along(uniq.aux.outlists))
      uniq.aux.outlists[[i]]$inputs[[4]] <- i-1

    o$model.aux["terms"] <- list(uniq.aux.outlists)

    # Remap client term pointers:
    for(i in seq_along(o$terms)){
      trm <- o$terms[[i]]
      if(is.null(trm$auxiliaries) || (naux <- length(list_rhs.formula(trm$auxiliaries)))==0) next
      trm$inputs[3+seq_len(naux)] <- omap[trm$inputs[3+seq_len(naux)]+1]
      o$terms[[i]] <- trm
    }
    o$slots.extra.aux <-  map(o$slots.extra.aux, ~omap[.+1])
    
    for(i in seq_along(m$terms)){
      trm <- m$terms[[i]]
      if(is.null(trm$auxiliaries) || (naux <- length(list_rhs.formula(trm$auxiliaries)))==0) next
      trm$inputs[3+seq_len(naux)] <- mmap[trm$inputs[3+seq_len(naux)]+1]
      m$terms[[i]] <- trm
    }
    m$slots.extra.aux <- map(m$slots.extra.aux, ~mmap[.+1])
    
    for(name in c("coef.names",
                  "inputs",
                  "minval",
                  "maxval",
                  "terms",
                  "offset",
                  "term.skipped",
                  "slots.extra.aux"))
      o[[name]] <- c(o[[name]], m[[name]])
    o$duration <- max(o$duration, m$duration)
    o$formula <- append_rhs.formula(o$formula, m$formula, keep.onesided=TRUE)
  }

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
