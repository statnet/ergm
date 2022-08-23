#  File R/param_names.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
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
    if(is.na(offset)) rep(TRUE, length(tocount))
    else if(offset) tocount
    else if(!offset) !tocount

  if(canonical) unlist(lapply(object$terms, function(term) term$coef.names))[tocount]
  else unlist(lapply(object$terms, function(term) NVL(names(term$params),term$coef.names)))[tocount]
}
