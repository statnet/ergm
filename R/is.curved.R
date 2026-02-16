#  File R/is.curved.R in package ergm, part of the Statnet suite of packages
#  for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################
###################################################################
## This file has utilities whose primary purpose is examining or ##
## manipulating ERGM formulas.                                   ##
###################################################################



#' Testing for curved exponential family
#' 
#' These functions test whether an ERGM fit, formula, or some other
#' object represents a curved exponential family.
#' 
#' Curvature is checked by testing if all model parameters are canonical.
#' 
#' @param object An [`ergm`] object or an ERGM formula.
#' @param \dots Arguments passed on to lower-level functions.
#' @return \code{TRUE} if the object represents a
#' curved exponential family; \code{FALSE} otherwise.
#' @keywords model
#' @export 
is.curved<-function(object,...) UseMethod("is.curved")

#' @rdname is.curved
#' @description The method for `NULL` always returns `FALSE` by
#'   convention.
#' @export
is.curved.NULL <- function(object, ...) FALSE # By convention.

#' @describeIn ergm_model Tests whether the model is curved.
#' @export
is.curved.ergm_model <- function(object, ...){
  NVL3(object$etamap$curved, length(.) > 0,
       map(object$terms, "map") %>% map_lgl(is.null) %>% all() %>% `!`)
}

#' @rdname is.curved 
#' @template response
#' @param basis See [ergm()].
#' @export
is.curved.formula<-function(object,response=NULL,basis=NULL,...){
  # If basis is not null, replace network in formula by basis.
  # In either case, let nw be network object from formula.
  if(is.null(nw <- basis)) {
    nw <- ergm.getnetwork(object)
  }
  
  nw <- ensure_network(nw)
  ergm_preprocess_response(nw, response)
  
  m<-ergm_model(object, nw, ...)
  is.curved(m)
}

#' @rdname is.curved 
#' @export
is.curved.ergm<-function(object,...){
  length(object$etamap$curved)>0
}
