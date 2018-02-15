#' Length of the parameter vector associated with an object or with its terms.
#'
#' This is a generic that returns the number of parameters associated with a model or a model fit. 
#' 
#' @param object An object for which number of parameters is defined.
#' 
#' @export
nparam <- function(object, ...){
  UseMethod("nparam")
}

#' @describeIn nparam
#'
#' By default, the length of the [coef()] vector is returned.
nparam.default <- function(object, ...){
  length(coef(object))
}

#' @describeIn nparam
#' A method to return the number of coefficients of an initialized [`ergm_model`].
#' 
#' @param offset If `NA` (the default), all model terms are counted;
#'   if \code{TRUE}, only offset terms are counted; and if
#'   \code{FALSE}, offset terms are skipped.
#' @param byterm Whether to return a vector of the numbers of
#'   coefficients for each term.
#' @template canonical
#' @export
nparam.ergm_model <- function(object, offset=NA, byterm=FALSE, canonical=FALSE, ...){
  terms <-
    if(is.na(offset)) object$terms
    else if(offset) object$terms[object$etamap$offset]
    else if(!offset) object$terms[!object$etamap$offset]
  
  out <-
    if(canonical){
      sapply(terms, function(term){
        length(term$coef.names)
      })
    }else{
      sapply(terms, function(term){
        ## curved term
        if(!is.null(term$params)) length(term$params)
        ## linear term
        else length(term$coef.names)
      })
    }
  if(byterm) out else sum(out)
}

#' @describeIn nparam
#' A method to return the number of parameters of an [`ergm`] fit.
#' 
#' @export
nparam.ergm <- function(object, offset=NA, ...){
  if(is.na(offset)) length(object$etamap$offsettheta)
  else if(offset) sum(object$etamap$offsettheta)
  else if(!offset) sum(!object$etamap$offsettheta)
}
