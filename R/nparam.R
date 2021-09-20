#  File R/nparam.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################
#' Length of the parameter vector associated with an object or with its terms.
#'
#' This is a generic that returns the number of parameters associated with a model or a model fit. 
#' 
#' @param object An object for which number of parameters is defined.
#' @param ... Additional arguments to methods.
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

#' @describeIn ergm_model Number of parameters of the model.
#' 
#' @template offset..filter
#' @param byterm Whether to return a vector of the numbers of
#'   coefficients for each term.
# #' @template canonical // Documented in one of the other methods of ergm_model.
#' @export
nparam.ergm_model <- function(object, canonical=FALSE, offset=NA, byterm=FALSE, ...){
  tocount <- if(canonical) object$etamap$offsetmap else object$etamap$offsettheta
  tocount <-
    if(is.na(offset)) rep(TRUE, length(tocount))
    else if(offset) tocount
    else if(!offset) !tocount

  if(byterm){
    terms <- object$terms
    lens <-
      if(canonical) terms %>% map("coef.names") %>% lengths
      else terms %>% map_int(function(term){
        ## curved term
        if(!is.null(term$params)) length(term$params)
        ## linear term
        else length(term$coef.names)
      })

    tocount %>% split(factor(rep(seq_along(terms), lens),levels=seq_along(terms))) %>% map_int(sum)
  }else sum(tocount)
}

#' @describeIn nparam
#' A method to return the number of parameters of an [`ergm`] fit.
#' 
#' @template offset..filter
#'
#' @export
nparam.ergm <- function(object, offset=NA, ...){
  if(is.na(offset)) length(object$etamap$offsettheta)
  else if(offset) sum(object$etamap$offsettheta)
  else if(!offset) sum(!object$etamap$offsettheta)
}
