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

#' @describeIn param_names
#' A method to return the curved or canonical parameter names of an initialized [`ergm_model`].
#'
#' @template canonical
#' @export
param_names.ergm_model <- function(object, canonical=FALSE, ...){
    if(canonical) object$coef.names
    else unlist(lapply(object$terms, function(term) NVL(names(term$params),term$coef.names)))
}
