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
#' 
#' @param formula a formula of the form \code{network ~ model.term(s)}
#' @param term.options a list of optional settings such as calculation tuning options to be passed to the `InitErgmTerm` functions.
#' @param extra.aux a list of auxiliary request formulas required elsewhere.
#' @param object See specific method documentation.
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
#' Main method for constructing a model object from an [ergm()]
#' formula of the form \code{network ~ model.term(s)} or \code{~
#' model.term(s)} with the network passed separately. Each term is
#' initialized via the \code{InitErgmTerm} functions to create a
#' \code{ergm_model} object.
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
ergm_model.formula <- function(object, nw, response=NULL, silent=FALSE, role="static",...,term.options=list(),extra.aux=list()){
  if (!is(object, "formula"))
    stop("Invalid model formula of class ",sQuote(class(object)),".", call.=FALSE)

  #' @importFrom statnet.common list_rhs.formula
  v<-list_rhs.formula(object)
  
  formula.env<-environment(object)
  
  model <- structure(list(formula=object, coef.names = NULL,
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
    # Now it is necessary to add the output to the model object
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
  
  termFun<-locate.InitFunction(term, paste0(termroot,"Term"), "ERGM term")  # check in all namespaces for function found anywhere

  # A kludge so that a term can read its own name.
  term.env <- new.env(parent=env) # Create a temporary new environment within env.
  assign(attr(termFun,"fname"), termFun, pos=term.env) # Assign the function body there.
  term<-as.call(list(as.name(attr(termFun,"fname"))))
  
  term[[2]] <- nw
  names(term)[2] <-  "nw"
  term[[3]] <- args
  names(term)[3] <- "arglist"
  dotdotdot <- c(if(!is.null(response)) list(response=response), list(role=role), c(term.options,list(...)))
  for(j in seq_along(dotdotdot)) {
    if(is.null(dotdotdot[[j]])) next
    term[[3+j]] <- dotdotdot[[j]]
    names(term)[3+j] <- names(dotdotdot)[j]
  }
  #Call the InitErgm function in the environment where the formula was created
  # so that it will have access to any parameters of the ergm terms
  out <- eval(term,term.env)
  # If SO package name not specified explicitly, autodetect.
  if(!is.null(out) && is.null(out$pkgname)) out$pkgname <- attr(termFun,"pkgname")
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
    model$duration <- c(model$duration,
                      NVL(outlist$duration, FALSE))
    model$terms[[termnumber]] <- outlist
  }
  model
}

