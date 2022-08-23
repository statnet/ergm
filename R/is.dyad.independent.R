#  File R/is.dyad.independent.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
###################################################################
## This file has utilities whose primary purpose is examining or ##
## manipulating ERGM formulas.                                   ##
###################################################################



#' Testing for dyad-independence
#' 
#' These functions test whether an ERGM fit, a formula, or some other
#' object represents a dyad-independent model.
#' 
#' Dyad independence is determined by checking if all of the
#' constituent parts of the object (formula, ergm terms, constraints,
#' etc.) are flagged as dyad-independent.
#' 
#' @param object The object to be tested for dyadic independence.
#' @param \dots Unused at this time.
#' @return \code{TRUE} if the model implied by the object is
#'   dyad-independent; \code{FALSE} otherwise.
#' @keywords model
#' @export is.dyad.independent
is.dyad.independent<-function(object,...) UseMethod("is.dyad.independent")

#' @rdname is.dyad.independent
#' @description The method for `NULL` always returns `TRUE` by
#'   convention.
#' @export
is.dyad.independent.NULL <- function(object, ...) TRUE # By convention.

#' @describeIn ergm_model Tests whether the model is dyad-independent.
#' @export
is.dyad.independent.ergm_model <- function(object, ...){
  ! any(sapply(object$terms, function(term) is.null(term$dependence) || term$dependence))
}

#' @rdname is.dyad.independent
#' @template response
#' @param basis See \code{\link{ergm}}.
#' @export
is.dyad.independent.formula<-function(object,response=NULL,basis=NULL,...){
  # If basis is not null, replace network in formula by basis.
  # In either case, let nw be network object from formula.
  if(is.null(nw <- basis)) {
    nw <- ergm.getnetwork(object)
  }
  
  nw <- as.network(nw)
  if(!is.network(nw)){
      stop("A network object on the LHS of the formula or via",
           " the 'basis' argument must be given")
    }

  ergm_preprocess_response(nw, response)
  m<-ergm_model(object, nw, ...)
  is.dyad.independent(m)
}

#' @rdname is.dyad.independent
#' @param object.obs For the [`ergm_conlist`] method, the observed data constraint.
#' @export
is.dyad.independent.ergm_conlist <- function(object, object.obs=NULL, ...){
  dind <- TRUE
  for(con in object){
    if(con$dependence) dind <- FALSE
  }
  dind && NVL3(object.obs, is.dyad.independent(.), TRUE)
}

#' @rdname is.dyad.independent
#' @export
is.dyad.independent.ergm<-function(object,...){
  NVL(object$MPLE_is_MLE,
      with(object,
           is.dyad.independent(formula,basis=network,...)
           && is.dyad.independent(object$constrained, object$constrained.obs))
      )
}
