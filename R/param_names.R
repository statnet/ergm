#  File R/param_names.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
#' Names of the parameters associated with an object.
#'
#' This is a generic that returns a vector giving the names of the parameters associated with a model or a model fit. 
#' 
#' @param object An object for which parameter names are defined.
#' @param ... Additional arguments to methods.
#' 
#' @export
param_names <- function(object, ...){
  UseMethod("param_names")
}

#' @describeIn param_names
#'
#' By default, the names of the [coef()] vector is returned.
#' @export
param_names.default <- function(object, ...){
  names(coef(object))
}

#' @describeIn ergm_model Parameter names of the model.
#'
#' @template canonical
#' @export
param_names.ergm_model <- function(object, canonical=FALSE, offset=NA, ...){
  tocount <- if(canonical) object$etamap$offsetmap else object$etamap$offsettheta
  tocount <-
    if(is.na(offset)) TRUE
    else if(offset) tocount
    else if(!offset) !tocount

  if(canonical) unlist(lapply(object$terms, function(term) term$coef.names))[tocount]
  else unlist(lapply(object$terms, function(term) NVL(names(term$params),term$coef.names)))[tocount]
}

#' @describeIn param_names a method for modifying parameter names of an object.
#'
#' @param value Specification for the new parameter names.
#'
#' @export
`param_names<-` <- function(object, ..., value){
  UseMethod("param_names<-")
}

#' @describeIn ergm_model Rename the parameters.
#'
#' @param value For [param_names<-()], either a character vector equal
#'   in length to the number of parameters of the specified type
#'   (though recycled as needed), or a [`list`] of two character
#'   vectors, one for non-canonical, the other for canonical, in which
#'   case `canonical=` will be ignored. `NA` elements preserve
#'   existing name.
#'
#' @export
`param_names<-.ergm_model` <- function(object, canonical = FALSE, ..., value){
  if(is.list(value)){
    if(length(value) != 2) stop("a character vector of a list of two character vectors is expected")
    param_names(object, FALSE) <- value[[1]]
    param_names(object, TRUE) <- value[[2]]
    return(object)
  }
  lens <- nparam(object, canonical, byterm = TRUE)
  value <- rep_len(value, sum(lens))

  values <- split(value, factor(rep(seq_along(lens), lens), levels = seq_along(lens)))
  for(i in seq_along(object$terms)){
    if(lens[i]){
      if(canonical || is.null(object$terms[[i]]$params)) object$terms[[i]]$coef.names %<>% replace(!is.na(values[[i]]), na.omit(values[[i]]))
      else names(object$terms[[i]]$params) %<>% replace(!is.na(values[[i]]), na.omit(values[[i]]))
    }
  }

  object
}
