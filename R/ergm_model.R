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

#' Internal Functions for Querying, Validating and Extracting from ERGM
#' Formulas
#' 
#' These are all functions that are generally not called directly by
#' users, but may be employed by other depending packages.  The
#' \code{ergm_model} function parses the given formula, and
#' initiliazes each ergm term via the \code{InitErgmTerm} functions to
#' create a \code{ergm_model} object.
#' 
#' @param formula a formula of the form \code{network ~ model.term(s)}
#' @param form same as formula, a formula of the form \code{network ~ model.term(s)}
#' @param object formula object to be updated
#' @param new new formula to be used in updating
#' @param from.new logical or character vector of variable names. controls how
#' environment of formula gets updated.
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
#' @note These functions are not meant to be called by the end-user but may be useful to extension developers. This API is not to be considered fixed and may change between versions.
#' @export
ergm_model <- function(object, ...){
  UseMethod("ergm_model")
}
#' @describeIn ergm_model
#'
#' Main method for constructing a model object from an [ergm()] formula of the form \code{network ~ model.term(s)}.
#' 
#' @param nw the network of interest
#' @param response charcter, name of edge attribute containing edge weights
#' @param silent logical, whether to print the warning messages from the
#' initialization of each model term; default=FALSE
#' @param role A hint about how the model will be used. Used primarily for
#' dynamic network models.
#' @param \dots additional parameters for model formulation
#'
#' @export
ergm_model.formula <- function(formula, nw, response=NULL, silent=FALSE, role="static",...) {
  if (!is(formula, "formula"))
    stop("Invalid model formula of class ",sQuote(class(formula)),".", call.=FALSE)

  if (length(formula) < 3) 
    stop("Model formula must have a left-hand-side.", call.=FALSE)

  #' @importFrom statnet.common list_rhs.formula
  v<-list_rhs.formula(formula)
  
  formula.env<-environment(formula)
  
  model <- structure(list(formula=formula, coef.names = NULL,
                      offset = NULL,
                      terms = NULL, networkstats.0 = NULL, etamap = NULL),
                 class = "ergm_model")

  termroot<-NVL2(response, "InitWtErgm", "InitErgm")

  
  for (i in 1:length(v)) {
    term <- v[[i]]
    
    if (is.call(term) && term[[1]] == "offset"){ # Offset term
      term <- term[[2]]
      model$offset <- c(model$offset,TRUE)
    }else{
      model$offset <- c(model$offset,FALSE)
    }
    ## term is now a call or a name that is not "offset".
    
    if(is.call(term)) { # This term has some arguments; save them.
      args <- term
      args[[1]] <- as.name("list")
    }else args <- list()
    
    termFun<-locate.InitFunction(term, paste0(termroot,"Term"), "ERGM term")  # check in all namespaces for function found anywhere
    
    term<-as.call(list(termFun))
    
    term[[2]] <- nw
    names(term)[2] <-  ""
    term[[3]] <- args
    names(term)[3] <- ""
    dotdotdot <- c(if(!is.null(response)) list(response=response), list(role=role), list(...))
    for(j in seq_along(dotdotdot)) {
      if(is.null(dotdotdot[[j]])) next
      term[[3+j]] <- dotdotdot[[j]]
      names(term)[3+j] <- names(dotdotdot)[j]
    }
    #Call the InitErgm function in the environment where the formula was created
    # so that it will have access to any parameters of the ergm terms
    outlist <- eval(term,formula.env)
    # If initialization fails without error (e.g., all statistics have been dropped), continue.
    if(is.null(outlist)){
      if(!silent) message("Note: Term ", deparse(v[[i]])," skipped because it contributes no statistics.")
      model$term.skipped <- c(model$term.skipped, TRUE)
      next
    }else model$term.skipped <- c(model$term.skipped, FALSE)
    # If SO package name not specified explicitly, autodetect.
    if(is.null(outlist$pkgname)) outlist$pkgname <- environmentName(environment(termFun))
    # If the term is an offset, rename the coefficient names and parameter names
    if(model$offset[length(model$offset)]){
      outlist$coef.names <- paste0("offset(",outlist$coef.names,")")
      if(!is.null(outlist$params))
        names(outlist$params) <- paste0("offset(",outlist$params,")")
    }
    # Now it is necessary to add the output to the model object
    model <- updatemodel.ErgmTerm(model, outlist)
  } 
  model$etamap <- ergm.etamap(model)

  # I.e., construct a vector of package names associated with the model terms.
  # Note that soname is not the same, since it's not guaranteed to be a loadable package.
  ergm.MCMC.packagenames(unlist(sapply(model$terms, "[[", "pkgname")))
  
  class(model) <- "ergm_model"
  model
}

#' @rdname ergm_model
#'
#' @description `ergm_model` is a deprecated name for the `ergm_model.formula` method.
#' @export ergm_model
ergm_model <- function(object, ...){
  .Deprecated("ergm_model")
  UseMethod("ergm_model")
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
    outlist$inputs <- c(ifelse(is.null(tmp), 0, tmp),
                        length(outlist$coef.names), 
                        length(outlist$inputs), outlist$inputs)
    model$minval <- c(model$minval,
                      rep(NVL(outlist$minval, -Inf),
                          length.out=length(outlist$coef.names)))
    model$maxval <- c(model$maxval,
                      rep(NVL(outlist$maxval, +Inf),
                          length.out=length(outlist$coef.names)))
    model$duration <- c(model$duration,
                      NVL(outlist$duration, FALSE))
    model$terms[[termnumber]] <- outlist
  }
  model
}

